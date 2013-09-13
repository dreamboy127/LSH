#include <stdio.h>
#include <stdlib.h>
#include "src/p_stable_LSH.h"

#define LSHDATABASE "bin/LSHDataBase.lshdb"
#define FEATUREBASE "bin/FeatureBase.lshfb"

int main(int argc,char* argv[])
{
	int HashNum=10;//L
	int VectorLength=4;//k
	int HASHSIZE;
	int DataBaseFeatureNum=1000;
	int dimention=784;
	int FeatureNum=9;
	double **feature;
	double **QueryPoint;
	int i;
	int j;
	int aclk;
	int bclk;
	FILE* OUT_fp;

	//FILE* DATABASE_FP=fopen("test/Pitch.dat","rt");
	//FILE* QUERY_FP=fopen("test/queryPitch.dat","rt");
	FILE* DATABASE_FP=fopen("test/mnist1k.dts","rt");
	FILE* QUERY_FP=fopen("test/mnist1k.q","rt");
	int ret;
//database
	
	feature=(double**)malloc(DataBaseFeatureNum*sizeof(double*));
	for (i=0;i<DataBaseFeatureNum;i++)
	{	
		feature[i]=(double*)malloc(dimention*sizeof(double));
		for(j=0;j<dimention;j++)
			fscanf(DATABASE_FP,"%lf",&feature[i][j]);
	}
	
	fclose(DATABASE_FP);
	
	HASHSIZE=DataBaseFeatureNum;
	
	aclk=clock();
	ret=LSH_Creat(LSHDATABASE,FEATUREBASE,feature,HashNum,VectorLength,HASHSIZE,DataBaseFeatureNum,dimention);
	bclk=clock();
	printf("LSH_Creat time:%lf\n",(bclk-aclk)/1000.0);
	if (ret<0)
	{
		printf("Creat Database Failed!\n");
		return 0;
	}
	for (i=0;i<DataBaseFeatureNum;i++)
	{	
		free(feature[i]);
	}
	free(feature);
	
//
	ret=LSH_Load(LSHDATABASE,FEATUREBASE,&HashNum,&VectorLength,&HASHSIZE,&DataBaseFeatureNum,&dimention);
	if (ret<0)
	{
		printf("Load Database Failed!\n");
		return 0;
	}
	QueryPoint=(double**)malloc(FeatureNum*sizeof(double*));
	for (i=0;i<FeatureNum;i++)
	{	
		QueryPoint[i]=(double*)malloc(dimention*sizeof(double));
		for(j=0;j<dimention;j++)
			fscanf(QUERY_FP,"%lf",&QueryPoint[i][j]);
	}
	fclose(QUERY_FP);
	//OUT_fp=fopen("out.txt","wt");
	OUT_fp=stdout;
	aclk=clock();
	LSH_MultiQuery(QueryPoint,FeatureNum,HashNum,VectorLength,HASHSIZE,DataBaseFeatureNum,dimention,OUT_fp);
	bclk=clock();
	//fclose(OUT_fp);
	printf("search time:%lf\n",(bclk-aclk)/1000.0);

	
	for (i=0;i<FeatureNum;i++)
	{	
		free(QueryPoint[i]);
	}
	free(QueryPoint);
	LSH_Release(HashNum,VectorLength,HASHSIZE,DataBaseFeatureNum);
	return 0;
}