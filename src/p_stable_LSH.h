#ifndef P_STABLE_LSH_H
#define P_STABLE_LSH_H

#include <stdio.h>
#define  PI 3.1415926

#define R 1
#define CandidateNum 3

typedef struct Point
{
	int pointid;
	struct Point* next;
}Point;
typedef struct HashBucket
{
	int pointNum;
	Point* point;
}HashBucket,HashBucketLoad,HashBucketCreat;

typedef struct QueryResultBucket
{
	int PointId;
	double distance;
}QueryResultBucket;

HashBucketCreat **HashsCreate;
HashBucketLoad **HashsLoad;
double** DataBaseFeature;


void InitLSH(int hashNum,int VectorLength,int HASHSIZE,int dimention);
void unInitLSH(int hashNum,int VectorLength);
int LSH_Creat(char* LSHDataBase,char* FeatureBase,double **feature,int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum,int dimention);
int LSH_Load(char* LSHDataBase,char* FeatureBase,int* HashNum,int *VectorLength,int *HASHSIZE,int* DataBaseFeatureNum,int* dimention);
void LSH_Release(int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum);
int LSH_SingleQuery(double *feature,int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum,int dimention,QueryResultBucket* QueryHash);
int LSH_MultiQuery(double **feature,int FeatureNum,int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum,int dimention,FILE* out_fp);
#endif