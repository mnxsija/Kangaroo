#include "Kangaroo.h"
#include <fstream>
#include "SECPK1/IntGroup.h"
#include "Timer.h"
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#ifndef WIN64
#include <dirent.h>
#include <pthread.h>
#endif

using namespace std;

bool Kangaroo::MergeWork(std::string& file1, std::string& file2, std::string& dest, bool printStat) {

  if(IsDir(file1) && IsDir(file2)) {
    return MergeWorkPartPart(file1, file2);
  }

  if(IsDir(file1)) {
    return MergeWorkPart(file1, file2, true);
  }

  if(dest.length() == 0) {
    ::printf("MergeWork: destination argument missing\n");
    return true;
  }

  double t0 = Timer::get_tick();
  uint32_t v1, v2;

  // Read headers and global parameters
  FILE* f1 = ReadHeader(file1, &v1, HEADW);
  FILE* f2 = ReadHeader(file2, &v2, HEADW);
  if(f1 == NULL || f2 == NULL) {
    if(f1) fclose(f1);
    if(f2) fclose(f2);
    return false;
  }

  if(v1 != v2) {
    ::printf("MergeWork: cannot merge workfiles of different versions\n");
    fclose(f1);
    fclose(f2);
    return true;
  }

  // Read global parameters
  uint32_t dp1, dp2;
  Int RS1, RE1, RS2, RE2;
  uint64_t count1, count2;
  double time1, time2;
  int keyCount1, keyCount2;

  ::fread(&dp1, sizeof(uint32_t), 1, f1);
  ::fread(&RS1.bits64, 32, 1, f1); RS1.bits64[4] = 0;
  ::fread(&RE1.bits64, 32, 1, f1); RE1.bits64[4] = 0;
  ::fread(&keyCount1, sizeof(int), 1, f1);

  ::fread(&dp2, sizeof(uint32_t), 1, f2);
  ::fread(&RS2.bits64, 32, 1, f2); RS2.bits64[4] = 0;
  ::fread(&RE2.bits64, 32, 1, f2); RE2.bits64[4] = 0;
  ::fread(&keyCount2, sizeof(int), 1, f2);

  if(!RS1.IsEqual(&RS2) || !RE1.IsEqual(&RE2)) {
    ::printf("MergeWork: File ranges differ\n");
    fclose(f1);
    fclose(f2);
    return true;
  }

  if(keyCount1 != keyCount2) {
    ::printf("MergeWork: Files have different number of keys\n");
    fclose(f1);
    fclose(f2);
    return true;
  }

  vector<Point> keys1, keys2;
  for(int i = 0; i < keyCount1; i++) {
    Point k;
    ::fread(&k.x.bits64, 32, 1, f1); k.x.bits64[4] = 0;
    ::fread(&k.y.bits64, 32, 1, f1); k.y.bits64[4] = 0;
    k.z.SetInt32(1);
    if(!secp->EC(k)) {
      ::printf("MergeWork: key in file1 does not lie on elliptic curve\n");
      fclose(f1);
      fclose(f2);
      return true;
    }
    keys1.push_back(k);

    ::fread(&k.x.bits64, 32, 1, f2); k.x.bits64[4] = 0;
    ::fread(&k.y.bits64, 32, 1, f2); k.y.bits64[4] = 0;
    k.z.SetInt32(1);
    if(!secp->EC(k)) {
      ::printf("MergeWork: key in file2 does not lie on elliptic curve\n");
      fclose(f1);
      fclose(f2);
      return true;
    }
    keys2.push_back(k);

    if(!keys1[i].equals(keys2[i])) {
      ::printf("MergeWork: keys at index %d differ\n", i);
      fclose(f1);
      fclose(f2);
      return true;
    }
  }

  ::fread(&count1, sizeof(uint64_t), 1, f1);
  ::fread(&time1, sizeof(double), 1, f1);
  ::fread(&count2, sizeof(uint64_t), 1, f2);
  ::fread(&time2, sizeof(double), 1, f2);

  ::printf("File %s: [DP%d]\n", file1.c_str(), dp1);
  ::printf("File %s: [DP%d]\n", file2.c_str(), dp2);

  // Set starting parameters
  endOfSearch = false;
  keysToSearch = keys1;
  keyIdx = 0;
  collisionInSameHerd = 0;
  rangeStart.Set(&RS1);
  rangeEnd.Set(&RE1);
  InitRange();
  InitSearchKey();

  // Open output file
  string tmpName = dest + ".tmp";
  FILE* f = fopen(tmpName.c_str(), "wb");
  if(f == NULL) {
    ::printf("\nMergeWork: Cannot open %s for writing\n", tmpName.c_str());
    ::printf("%s\n", ::strerror(errno));
    fclose(f1);
    fclose(f2);
    return true;
  }

  dpSize = (dp1 < dp2) ? dp1 : dp2;
  if(!SaveHeader(tmpName, f, HEADW, count1 + count2, time1 + time2)) {
    fclose(f1);
    fclose(f2);
    fclose(f);
    return true;
  }

  ::printf("Merging");

  uint64_t nbDP = 0;
  uint32_t hDP;
  uint32_t hDuplicate;
  Int d1, d2;
  uint32_t type1, type2;

  for(uint32_t h = 0; h < HASH_SIZE && !endOfSearch; h++) {
    if(h % (HASH_SIZE / 64) == 0) ::printf(".");

    int mStatus = HashTable::MergeH(h, f1, f2, f, &hDP, &hDuplicate, &d1, &type1, &d2, &type2);
    switch(mStatus) {
      case ADD_OK:
        break;
      case ADD_COLLISION:
        CollisionCheck(&d1, type1, &d2, type2);
        break;
    }

    nbDP += hDP;
    collisionInSameHerd += hDuplicate;
  }

  fclose(f1);
  fclose(f2);
  fclose(f);

  double t1 = Timer::get_tick();

  if(!endOfSearch) {
    remove(dest.c_str());
    rename(tmpName.c_str(), dest.c_str());
    ::printf("Done [%s]\n", GetTimeStr(t1 - t0).c_str());
  } else {
    remove(tmpName.c_str());
    return true;
  }

  if(printStat) {
    ::printf("Dead kangaroo: %" PRId64 "\n", collisionInSameHerd);
    ::printf("Total f1+f2: DP count 2^%.2f\n", log2((double)nbDP));
  } else {
    offsetTime = time1 + time2;
    offsetCount = count1 + count2;
  }

  return false;
}

typedef struct File {
  std::string name;
  uint64_t size;
} File;

bool sortBySize(const File& lhs,const File& rhs) { return lhs.size > rhs.size; }

void Kangaroo::MergeDir(std::string& dirName,std::string& dest) {

  vector<File> listFiles;

#ifdef WIN64
  WIN32_FIND_DATA ffd;
  HANDLE hFind;

  hFind = FindFirstFile((dirName.c_str()+string("\\*")).c_str(),&ffd);
  if( hFind==INVALID_HANDLE_VALUE ) {
    ::printf("FindFirstFile Error: %d\n",GetLastError());
    return;
  }

  do {
    if((ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)==0) {
      uint32_t version;
      string fName = dirName.c_str() + string("\\") + string(ffd.cFileName);
      FILE *f = ReadHeader(fName,&version,HEADW);
      if(f) {
        File e;
        e.name = fName;
        _fseeki64(f,0,SEEK_END);
        e.size = (uint64_t)_ftelli64(f);
        listFiles.push_back(e);
        fclose(f);
      }
    }
  } while(FindNextFile(hFind,&ffd) != 0);
  FindClose(hFind);

#else

  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(dirName.c_str())) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      if ( ent->d_type != 0x8) continue;
      uint32_t version;
      string fName = dirName + "/" + string(ent->d_name);
      FILE *f = ReadHeader(fName,&version,HEADW);
      if(f) {
        File e;
        e.name = fName;
        fseeko(f,0,SEEK_END);
        e.size = (uint64_t)ftello(f);
        listFiles.push_back(e);
        fclose(f);
      }

    }
    closedir(dir);
  } else {
    ::printf("opendir(%s) Error:\n",dirName.c_str());
    perror("");
    return;
  }

#endif

  std::sort(listFiles.begin(),listFiles.end(),sortBySize);
  int lgth = (int)listFiles.size();

  if(IsDir(dest)==1) {

    // Partitioned merge
    bool end = false;
    for(int i = 0; i < lgth && !end; i++) {
      ::printf("\n## File #%d/%d\n",i+1,lgth);
      end = MergeWorkPart(dest,listFiles[i].name,i == lgth - 1);
    }

  } else {

    // Standard merge
    if(listFiles.size() < 2) {
      ::printf("MergeDir: less than 2 work files in the directory\n");
      return;
    }

    int i = 0;
    ::printf("\n## File #1/%d\n",lgth - 1);
    bool end = MergeWork(listFiles[0].name,listFiles[1].name,dest,lgth == 2);
    for(int i = 2; i < lgth && !end; i++) {
      ::printf("\n## File #%d/%d\n",i,lgth - 1);
      end = MergeWork(dest,listFiles[i].name,dest,i == lgth - 1);
    }

  }

 
}
