#include "p_stable_LSH.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

double r=4.0;
double **b;
double ***a;
int **weight1;
//int **weight2;

// 4294967291 = 2^32-5
#define UH_PRIME_DEFAULT 4294967291U
// 2^29
#define MAX_HASH_RND 536870912U
// 2^32-1
#define TWO_TO_32_MINUS_1 4294967295U

double EuclideanDistance(int dimension, double* p1, double* p2){
	double result = 0;
	int i;
	for (i = 0; i < dimension; i++){
		result += (p1[i] - p2[i])*(p1[i] - p2[i]);
	}
	
	return sqrt(result);
}
/************************************************************************/
/* qsort compare����                                                    */
/************************************************************************/
/*int compare_entries(const void *a, const void *b){
	return (((QueryHashBucket*)a)->entries > ((QueryHashBucket*)b)->entries)? 1:-1;
}*/

double AverageRandom(double min,double max)//ƽ���ֲ�
{
    int minInteger = (int)(min*10000);
    int maxInteger = (int)(max*10000);
    int randInteger = rand()*rand();
    int diffInteger = maxInteger - minInteger;
    int resultInteger = randInteger % diffInteger + minInteger;
    return resultInteger/10000.0;
}

double Normal(double x,double miu,double sigma) //�����ܶȺ���
{
	return 1.0/sqrt(2*PI*sigma) * exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}
double NormalRandom(double miu,double sigma,double min,double max)//������̬�ֲ������
{
    double x;
    double dScope;
    double y;
    do
	{
		x = AverageRandom(min,max); 
        y = Normal(x, miu, sigma);
        dScope = AverageRandom(0, Normal(miu,miu,sigma));
	}while( dScope > y);
	
	return x;
}
#define RAND_MAX_64 0x7fffffff
int Random32(int rangeStart, int rangeEnd)
{
	int r;
	if (RAND_MAX_64 >= rangeEnd - rangeStart) {
		r = rangeStart + (int)((rangeEnd - rangeStart + 1.0) * rand() / (RAND_MAX_64 + 1.0));
	} else {
		r = rangeStart + (int)((rangeEnd - rangeStart + 1.0) * ((LONG64)rand() * ((LONG64)RAND_MAX_64 + 1) + (LONG64)rand()) / ((LONG64)RAND_MAX_64 * ((LONG64)RAND_MAX_64 + 1) + (LONG64)RAND_MAX + 1.0));
	}
	return r;
}
/************************************************************************/
/*  ��ʼ��LSH��r,b,��ϣ��������г�ʼ�� r,b,a 
	����hashNum����ϣ�� ÿ�����СΪHASHSIZE
	����ά�ȴ�СΪdimention	
	���ڴ������ݿ����										            */
/************************************************************************/
void InitLSH(int hashNum,int VectorLength,int HASHSIZE,int dimention)
{
	int j,k,i;
	srand(time(NULL));
	
	r=4.0;//LSH�������r
	//b=AverageRandom(0,r);//LSH�������b
	a=(double***)malloc(sizeof(double**)*hashNum);
	b=(double**)malloc(sizeof(double*)*hashNum);
	weight1=(int**)malloc(sizeof(int*)*hashNum);
	//weight2=(int**)malloc(sizeof(int*)*hashNum);
	//a=(double**)malloc(sizeof(double*)*hashNum);
	HashsCreate=(HashBucketCreat**)malloc(sizeof(HashBucketCreat*)*hashNum);
	for (j=0;j<hashNum;j++)
	{
		a[j]=(double**)malloc(sizeof(double*)*VectorLength);
		b[j]=(double*)malloc(sizeof(double)*VectorLength);
		weight1[j]=(int*)malloc(sizeof(int)*VectorLength);
		//weight2[j]=(int*)malloc(sizeof(int)*VectorLength);
		for (k=0;k<VectorLength;k++)
		{
			a[j][k]=(double*)malloc(sizeof(double)*dimention);
			for (i=0;i<dimention;i++)
			{
				a[j][k][i]=NormalRandom(0,1,-1.5,1.5);
			}
			b[j][k]=AverageRandom(0,r);
			weight1[j][k]=Random32(1,MAX_HASH_RND);
			//weight2[j][k]=Random32(1,MAX_HASH_RND);
		}
		
		HashsCreate[j]=(HashBucketCreat*)malloc(sizeof(HashBucketCreat)*HASHSIZE);
		memset(HashsCreate[j],0,sizeof(HashBucketCreat)*HASHSIZE);
	}
}
/************************************************************************/
/*  LSH����һЩ�ڴ�����ͷ�		
	���ڴ������ݿ����										            */
/************************************************************************/
void unInitLSH(int hashNum,int VectorLength)
{
	int i,j;
	for (j=0;j<hashNum;j++)
	{
		for (i=0;i<VectorLength;i++)
		{
			free(a[j][i]);
		}
		free(a[j]);
		free(b[j]);
		free(weight1[j]);
		//free(weight2[j]);
		free(HashsCreate[j]);
	}
	free(a);
	free(b);
	free(weight1);
	//free(weight2);
	free(HashsCreate);
}
/************************************************************************/
/*  �����������������ϣֵ��ʹ���˹�ϣ����a,b,r				            */
/************************************************************************/
unsigned int HashFamilyFun(double *feature,double *a_temp,double b_temp,double r_temp,int dimention)//��ϣ����
{
	double result=b_temp;
	int i;
	for(i=0;i<dimention;i++)
	{
		result+=feature[i]*(*(a_temp+i));
	}
	return (unsigned int)floor(result/r_temp);//���ع�ϣ���
}
/************************************************************************/
/*  ������������Ӧ��ֵ���ݹ�ϣֵ���뵽��ϣ��Ķ�Ӧλ��
	i�����i����ϣ��
	hashval�����ϣֵ
	Pointid��������������Ӧ��ֵ��������Point��ID��
	���ڴ������ݿ����										            */
/************************************************************************/
int AddPoint2Hash(int i,int hashval1,int Pointid)
{
	Point** CurrentPoint=&HashsCreate[i][hashval1].point;
	HashsCreate[i][hashval1].pointNum++;
	while(*CurrentPoint!=NULL)
		CurrentPoint=&(*CurrentPoint)->next;
	*CurrentPoint=(Point*)malloc(sizeof(Point));
	if (*CurrentPoint==NULL)
	{
		return -1;
	}
	(*CurrentPoint)->pointid=Pointid;
	(*CurrentPoint)->next=NULL;
	return 0;
}
/************************************************************************/
/*  ��һ�������������뵽���й�ϣ����
	���ڴ������ݿ����										            */
/************************************************************************/
int AddPoint2Hashs(double *feature,int Pointid,int hashNum,int VectorLength,int HASHSIZE,int dimention)
{
	int i,j;
	LONG64 hashvalue1;
	//LONG64 hashvalue2;

	for (i=0;i<hashNum;i++)
	{
		hashvalue1=0;
		//hashvalue2=0;
		for (j=0;j<VectorLength;j++)
		{
			hashvalue1 += HashFamilyFun(feature,a[i][j],b[i][j],r,dimention)*weight1[i][j];
			hashvalue1 = (hashvalue1 & TWO_TO_32_MINUS_1) + 5 * (hashvalue1 >> 32);
			if (hashvalue1 >= UH_PRIME_DEFAULT) {
				hashvalue1 = hashvalue1 - UH_PRIME_DEFAULT;
			}
			/*hashvalue2 += HashFamilyFun(feature,a[i][j],b[i][j],r,dimention)*weight2[i][j];
			hashvalue2 = (hashvalue2 & TWO_TO_32_MINUS_1) + 5 * (hashvalue2 >> 32);
			if (hashvalue2 >= UH_PRIME_DEFAULT) {
				hashvalue2 = hashvalue2 - UH_PRIME_DEFAULT;
			}*/
		}
		
		hashvalue1=hashvalue1%HASHSIZE;
		if (AddPoint2Hash(i,(int)hashvalue1,Pointid)!=0)
		{
			return -1;
		}
	}
	return 0;
}
/************************************************************************/
/*  ��featurenum�������������뵽���й�ϣ����
	���ڴ������ݿ����										            */
/************************************************************************/
int LSH_Add_Feature(double **feature,int DataBaseFeatureNum,int HashNum,int VectorLength,int HASHSIZE,int dimention)
{
	int i;

	for (i=0;i<DataBaseFeatureNum;i++)
	{
		if (AddPoint2Hashs(feature[i],i+1,HashNum,VectorLength,HASHSIZE,dimention)!=0)
		{
			return -1;
		}
	}
	
	return 0;
}
/************************************************************************/
/*  �������ݿ���̣������������ݿ���һЩ����д�뵽LSHDataset�ļ���
	���ڴ������ݿ����										            */
/************************************************************************/
int LSH_Creat(char* LSHDataBase,char* FeatureBase,double **feature,int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum,int dimention)
{
	FILE* LSHDATABase_FP;
	FILE* FeatureBase_FP;
	int i,j;
	if (DataBaseFeatureNum<CandidateNum)
	{
		//���ݿ���������̫��
		return -1;
	}
	InitLSH(HashNum,VectorLength,HASHSIZE,dimention);

	if (LSH_Add_Feature(feature,DataBaseFeatureNum,HashNum,VectorLength,HASHSIZE,dimention)!=0)
	{
		return -1;
	}
	
	LSHDATABase_FP=fopen(LSHDataBase,"wb");

	fwrite(&HashNum,sizeof(int),1,LSHDATABase_FP);
	fwrite(&VectorLength,sizeof(int),1,LSHDATABase_FP);
	fwrite(&HASHSIZE,sizeof(int),1,LSHDATABase_FP);
	fwrite(&DataBaseFeatureNum,sizeof(int),1,LSHDATABase_FP);
	fwrite(&dimention,sizeof(int),1,LSHDATABase_FP);
	for (i=0;i<HashNum;i++)
	{
		for (j=0;j<VectorLength;j++)
		{
			fwrite(a[i][j],sizeof(double),dimention,LSHDATABase_FP);
			fwrite(&b[i][j],sizeof(double),1,LSHDATABase_FP);
			fwrite(&weight1[i][j],sizeof(int),1,LSHDATABase_FP);
			//fwrite(&weight2[i][j],sizeof(int),1,LSHDATABase_FP);
		}
	}
	fwrite(&r,sizeof(double),1,LSHDATABase_FP);
	for (i=0;i<HashNum;i++)
	{
		int j;
		for (j=0;j<HASHSIZE;j++)
		{
			int k;
			Point* CurrentPoint=HashsCreate[i][j].point;
			Point* tmpCurrentPoint;
			fwrite(&HashsCreate[i][j].pointNum,sizeof(int),1,LSHDATABase_FP);
			
			for (k=0;k<HashsCreate[i][j].pointNum;k++)
			{
				fwrite(&(CurrentPoint->pointid),sizeof(int),1,LSHDATABase_FP);
				tmpCurrentPoint=CurrentPoint;
				CurrentPoint=CurrentPoint->next;
				free(tmpCurrentPoint);
			}
		}
	}
	fclose(LSHDATABase_FP);
	FeatureBase_FP=fopen(FeatureBase,"wb");
	for (i=0;i<DataBaseFeatureNum;i++)
	{
		int j;
		for (j=0;j<dimention;j++)
		{
			fwrite(&feature[i][j],sizeof(double),1,FeatureBase_FP);
		}
	}
	fclose(FeatureBase_FP);
	unInitLSH(HashNum,VectorLength);
	return 0;
}
/************************************************************************/
/*  �������ݿ��ļ�����һЩ�������ϣ������
	����Ĺ�ϣ��ṹHashsLoad�ʹ������ݿ��е�HashsCreat�ṹ��ͬ
	����˳�ʼ�����������Խ��в�ѯ
	���ڲ�ѯ����											            */
/************************************************************************/
int LSH_Load(char* LSHDataBase,char* FeatureBase,int* HashNum,int *VectorLength,int *HASHSIZE,int* DataBaseFeatureNum,int* dimention)
{
	FILE* LSHDATASET_FP;
	FILE* FeatureBase_FP;
	int i,j;
	LSHDATASET_FP=fopen(LSHDataBase,"rb");
	fread(HashNum,sizeof(int),1,LSHDATASET_FP);
	fread(VectorLength,sizeof(int),1,LSHDATASET_FP);
	fread(HASHSIZE,sizeof(int),1,LSHDATASET_FP);
	fread(DataBaseFeatureNum,sizeof(int),1,LSHDATASET_FP);
	if (*DataBaseFeatureNum<CandidateNum)
	{
		//���ݿ���������̫��
		fclose(LSHDATASET_FP);
		return -1;
	}
	fread(dimention,sizeof(int),1,LSHDATASET_FP);
	a=(double***)malloc(sizeof(double**)**HashNum);
	b=(double**)malloc(sizeof(double*)**HashNum);
	weight1=(int**)malloc(sizeof(int*)**HashNum);
	//weight2=(int**)malloc(sizeof(int*)**HashNum);
	for (i=0;i<*HashNum;i++)
	{
		a[i]=(double**)malloc(sizeof(double*)**VectorLength);
		b[i]=(double*)malloc(sizeof(double)**VectorLength);
		weight1[i]=(int*)malloc(sizeof(int)**VectorLength);
		//weight2[i]=(int*)malloc(sizeof(int)**VectorLength);
		for(j=0;j<*VectorLength;j++)
		{
			a[i][j]=(double*)malloc(sizeof(double)**dimention);
			fread(a[i][j],sizeof(double),*dimention,LSHDATASET_FP);
			fread(&b[i][j],sizeof(double),1,LSHDATASET_FP);
			fread(&weight1[i][j],sizeof(int),1,LSHDATASET_FP);
			//fread(&weight2[i][j],sizeof(int),1,LSHDATASET_FP);
		}
	}
	fread(&r,sizeof(double),1,LSHDATASET_FP);
	HashsLoad=(HashBucketLoad**)malloc(sizeof(HashBucketLoad*)**HashNum);
	for (i=0;i<*HashNum;i++)
	{
		int j;
		HashsLoad[i]=(HashBucketLoad*)malloc(sizeof(HashBucketLoad)**HASHSIZE);
		for (j=0;j<*HASHSIZE;j++)
		{
			int k;
			Point** CurrentPoint;
			fread(&HashsLoad[i][j].pointNum,sizeof(int),1,LSHDATASET_FP);

			CurrentPoint=&HashsLoad[i][j].point;
			for (k=0;k<HashsLoad[i][j].pointNum;k++)
			{
				int tmppoint;
				fread(&tmppoint,sizeof(int),1,LSHDATASET_FP);

				*CurrentPoint=(Point*)malloc(sizeof(Point));
				if (*CurrentPoint==NULL)
				{
					return -1;
				}
				(*CurrentPoint)->pointid=tmppoint;
				(*CurrentPoint)->next=NULL;
				CurrentPoint=&(*CurrentPoint)->next;
			}
		}
	}
	
	fclose(LSHDATASET_FP);
	DataBaseFeature=(double**)malloc(sizeof(double*)**DataBaseFeatureNum);

	FeatureBase_FP=fopen(FeatureBase,"rb");
	for (i=0;i<*DataBaseFeatureNum;i++)
	{
		int j;
		DataBaseFeature[i]=(double*)malloc(sizeof(double)**dimention);
		for (j=0;j<*dimention;j++)
		{
			fread(&DataBaseFeature[i][j],sizeof(double),1,FeatureBase_FP);
		}
	}
	fclose(FeatureBase_FP);
	return 0;
}
/************************************************************************/
/*  �ͷ�LSH�е�һЩ�ڴ�
	���ڲ�ѯ����											            */
/************************************************************************/
void LSH_Release(int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum)
{
	int j,i;
	for (j=0;j<HashNum;j++)
	{
		for(i=0;i<VectorLength;i++)
		{
			free(a[j][i]);
		}
		free(a[j]);
		free(b[j]);
		free(weight1[j]);
		//free(weight2[j]);
		for (i=0;i<HASHSIZE;i++)
		{
			Point* CurrentPoint=HashsLoad[j][i].point;
			Point* TmpCurrentPoint;
			int k;
			for (k=0;k<HashsLoad[j][i].pointNum;k++)
			{
				TmpCurrentPoint=CurrentPoint;
				CurrentPoint=CurrentPoint->next;
				free(TmpCurrentPoint);
			}
		}
		free(HashsLoad[j]);
	}
	free(a);
	free(b);
	free(weight1);
	//free(weight2);
	free(HashsLoad);
	for (j=0;j<DataBaseFeatureNum;j++)
	{
		free(DataBaseFeature[j]);
	}
	free(DataBaseFeature);
}
/************************************************************************/
/*  ���ڵ������������Ĳ�ѯ����ѯ������뵽LSH_Query_Result��
	��ѯ���ΪǰCandidateNum��������С�Ĵ����QueryHash��
	���ڲ�ѯ����											            */
/************************************************************************/
int LSH_SingleQuery(double *feature,int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum,int dimention,QueryResultBucket* QueryResult)
{
	int i,j,k;
	LONG64 key1;
	//LONG64 key2;
	int* MarkPointTable=(int*)malloc(sizeof(int)*DataBaseFeatureNum);
	
	memset(MarkPointTable,0,sizeof(int)*DataBaseFeatureNum);
	for (i=0;i<CandidateNum;i++)
	{
		QueryResult[i].distance=R;
	}
	for (i=0;i<HashNum;i++)
	{
		Point* CurrentPoint;

		key1=0;
		//key2=0;
		for (j=0;j<VectorLength;j++)
		{
			key1 += HashFamilyFun(feature,a[i][j],b[i][j],r,dimention)*weight1[i][j];
			key1 = (key1 & TWO_TO_32_MINUS_1) + 5 * (key1 >> 32);
			if (key1 >= UH_PRIME_DEFAULT) {
				key1 = key1 - UH_PRIME_DEFAULT;
			}
			/*key2 += HashFamilyFun(feature,a[i][j],b[i][j],r,dimention)*weight2[i][j];
			key2 = (key2 & TWO_TO_32_MINUS_1) + 5 * (key2 >> 32);
			if (key2 >= UH_PRIME_DEFAULT) {
				key2 = key2 - UH_PRIME_DEFAULT;
			}*/
		}

		key1=key1%HASHSIZE;

		CurrentPoint=HashsLoad[i][key1].point;
		for (j=0;j<HashsLoad[i][key1].pointNum;j++)
		{
			int index=CurrentPoint->pointid;
			double Distance;
			if(MarkPointTable[index-1]==1)
			{
				continue;
			}
			MarkPointTable[index-1]=1;
			Distance=EuclideanDistance(dimention,DataBaseFeature[index-1],feature);
			for (k=0;k<CandidateNum;k++)
			{
				if (Distance>R||Distance>QueryResult[CandidateNum-1].distance)
				{
					break;
				}
				if (QueryResult[k].distance>Distance)
				{
					int n;
					for (n=CandidateNum-1;n>k;n--)
					{
						QueryResult[n].distance=QueryResult[n-1].distance;
						QueryResult[n].PointId=QueryResult[n-1].PointId;
					}
					QueryResult[k].distance=Distance;
					QueryResult[k].PointId=index;
					break;
				}
			}
			CurrentPoint=CurrentPoint->next;
		}
	}

	free(MarkPointTable);
	return 0;
}
/************************************************************************/
/*  ���ڶ������������
	���ڲ�ѯ����											            */
/************************************************************************/
int LSH_MultiQuery(double **feature,int FeatureNum,int HashNum,int VectorLength,int HASHSIZE,int DataBaseFeatureNum,int dimention,FILE* out_fp)
{
	int i,j;
	QueryResultBucket QueryResult[CandidateNum];
	
	for (i=0;i<FeatureNum;i++)
	{
		for (j=0;j<CandidateNum;j++)
		{
			QueryResult[j].PointId=0;
		}
		LSH_SingleQuery(feature[i],HashNum,VectorLength,HASHSIZE,DataBaseFeatureNum,dimention,QueryResult);
		
		fprintf(out_fp,"the %dth point Query Result:\n",i+1);
		for (j=0;j<CandidateNum;j++)
		{
			if (QueryResult[j].PointId==0)
			{
				break;
			}
			fprintf(out_fp,"%d:%lf\n",QueryResult[j].PointId,QueryResult[j].distance);
		}
	}
	
	return 0;
}