#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include "mpi.h"


using namespace std;

int prank;  
int pnum;  

const double infper = 50.0;
const double RandomMultiplier = 10;

int Min(int A, int B) {
	int Result = (A < B) ? A : B;
	if ((A < 0) && (B >= 0)) Result = B;
	if ((B < 0) && (A >= 0)) Result = A;
	if ((A < 0) && (B < 0)) Result = -1;
	return Result;
}

void copyM(int *pMatrix, int Size, int *pMatrixCopy) {
	copy(pMatrix, pMatrix + Size * Size, pMatrixCopy);
}
 
bool compareM(int *pMatrix1, int *pMatrix2, int Size) {
	return equal(pMatrix1, pMatrix1 + Size * Size, pMatrix2);
}
 
void SerialFloyd(int *pMatrix, int Size) {
	int t1, t2;
	for (int k = 0; k < Size; k++)
		for (int i = 0; i < Size; i++)
			for (int j = 0; j < Size; j++)
				if ((pMatrix[i * Size + k] != -1) &&
					(pMatrix[k * Size + j] != -1)) {
					t1 = pMatrix[i * Size + j];
					t2 = pMatrix[i * Size + k] + pMatrix[k * Size + j];
					pMatrix[i * Size + j] = Min(t1, t2);
				}
}
 
void printM(int *pMatrix, int RowCount, int ColCount) {

	for (int i = 0; i < RowCount; i++) {
		for (int j = 0; j < ColCount; j++) {
			printf("%7d", pMatrix[i * ColCount + j]);
			fflush(stdout);
		}
		printf("\n");
		fflush(stdout);
	}
}

void term(int *pMatrix, int *pProcRows) {
	if (prank == 0)
		delete[]pMatrix;
	delete[]pProcRows;
}
 
void initialezeR(int *pMatrix, int Size) {
	srand((unsigned)time(0));
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
			if (i != j) {
				if ((rand() % 100) < infper)
					pMatrix[i * Size + j] = -1;
				else
					pMatrix[i * Size + j] = rand() + 1;
			}
			else
				pMatrix[i * Size + j] = 0;
}

void dataDist(int *pMatrix, int *pProcRows, int Size, int RowNum) {
	int *pSendNum;  
	int *pSendInd;  
	int RestRows = Size;  
	pSendInd = new int[pnum];
	pSendNum = new int[pnum];
	RowNum = Size / pnum;
	pSendNum[0] = RowNum * Size;
	pSendInd[0] = 0;
	for (int i = 1; i < pnum; i++) {
		RestRows -= RowNum;
		RowNum = RestRows / (pnum - i);
		pSendNum[i] = RowNum * Size;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}
	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_INT,
		pProcRows, pSendNum[prank], MPI_INT, 0, MPI_COMM_WORLD);
	delete[]pSendNum;
	delete[]pSendInd;
}
 
void ResCollect(int *pMatrix, int *pProcRows, int Size, int RowNum) {
	int *pReceiveNum;  
	int *pReceiveInd;  
	int RestRows = Size;  
	pReceiveNum = new int[pnum];
	pReceiveInd = new int[pnum];
	RowNum = Size / pnum;
	pReceiveInd[0] = 0;
	pReceiveNum[0] = RowNum * Size;
	for (int i = 1; i < pnum; i++) {
		RestRows -= RowNum;
		RowNum = RestRows / (pnum - i);
		pReceiveNum[i] = RowNum * Size;
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
	 
	MPI_Gatherv(pProcRows, pReceiveNum[prank], MPI_INT,
		pMatrix, pReceiveNum, pReceiveInd, MPI_INT, 0, MPI_COMM_WORLD);
	delete[]pReceiveNum;
	delete[]pReceiveInd;
}


 
void RowDist(int *pProcRows, int Size, int RowNum, int k, int
	*pRow) {
	int ProcRowRank;  
	int ProcRowNum;  
	int RestRows = Size;
	int Ind = 0;
	int Num = Size / pnum;
	for (ProcRowRank = 1; ProcRowRank < pnum + 1; ProcRowRank++) {
		if (k < Ind + Num) break;
		RestRows -= Num;
		Ind += Num;
		Num = RestRows / (pnum - ProcRowRank);
	}
	ProcRowRank = ProcRowRank - 1;
	ProcRowNum = k - Ind;
	if (ProcRowRank == prank) 
		copy(&pProcRows[ProcRowNum*Size], &pProcRows[(ProcRowNum + 1)*Size], pRow);
	MPI_Bcast(pRow, Size, MPI_INT, ProcRowRank, MPI_COMM_WORLD);
}

void ParallelFloyd(int *pProcRows, int Size, int RowNum) {
	int *pRow = new int[Size];
	int t1, t2;
	for (int k = 0; k < Size; k++) { 
		RowDist(pProcRows, Size, RowNum, k, pRow);
		for (int i = 0; i < RowNum; i++)
			for (int j = 0; j < Size; j++)
				if ((pProcRows[i * Size + k] != -1) &&
					(pRow[j] != -1)) {
					t1 = pProcRows[i * Size + j];
					t2 = pProcRows[i * Size + k] + pRow[j];
					pProcRows[i * Size + j] = Min(t1, t2);
				}
	}
	delete[]pRow;
}

void printMpar(int *pProcRows, int Size, int RowNum) {
	for (int i = 0; i < pnum; i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (prank == i) {
			printf("Process Rank = %d\n", prank);
			fflush(stdout);
			printf("Process rows:\n");
			fflush(stdout);
			printM(pProcRows, RowNum, Size);
			printf("\n");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
 
void testDist(int *pMatrix, int *pProcRows, int Size, int RowNum) {
	MPI_Barrier(MPI_COMM_WORLD);
	if (prank == 0) {
		printf("Initial adjacency matrix:\n");
		printM(pMatrix, Size, Size);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	printMpar(pProcRows, Size, RowNum);
}
 
void initializeP(int *&pMatrix, int *&pProcRows, int& Size, int&
	RowNum) {
	setvbuf(stdout, 0, _IONBF, 0);
	if (prank == 0) {
		do {
			printf("Enter the number of vertices: ");
			scanf_s("%d", &Size);
			if (Size < pnum)
				printf("The number of vertices should be greater then number of processes\n");
		} while (Size < pnum);	 
	}
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int RestRows = Size;
	for (int i = 0; i < prank; i++)
		RestRows = RestRows - RestRows / (pnum - i);
	RowNum = RestRows / (pnum - prank);
	pProcRows = new int[Size * RowNum];
	if (prank == 0) {
		pMatrix = new int[Size * Size];
		initialezeR(pMatrix, Size);
	}
}


int main(int argc, char* argv[]) {
	int *pMatrix;  
	int Size;  
	int *pProcRows;  
	int RowNum;  
	double start, finish;
	double duration = 0.0;
	int *pSerialMatrix = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &prank);
	if (prank == 0)
		printf("Parallel Floyd algorithm\n");
	initializeP(pMatrix, pProcRows, Size, RowNum);
	if (prank == 0) { 
		pSerialMatrix = new int[Size * Size];
		copyM(pMatrix, Size, pSerialMatrix);
	}
	start = MPI_Wtime();
	dataDist(pMatrix, pProcRows, Size, RowNum);
	testDist(pMatrix, pProcRows, Size, RowNum);
	ParallelFloyd(pProcRows, Size, RowNum);
	printMpar(pProcRows, Size, RowNum); 
	ResCollect(pMatrix, pProcRows, Size, RowNum);
	if (prank == 0)
		printM(pMatrix, Size, Size);
	finish = MPI_Wtime();
	duration = finish - start;
	if (prank == 0)
		printf("Time of execution: %f\n", duration);
	if (prank == 0)
		delete[]pSerialMatrix;
	term(pMatrix, pProcRows);
	MPI_Finalize();

	system("pause");

	return 0;
}



