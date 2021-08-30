#include <cstdlib>
#include <cstdio>
#include <ctime>

using namespace std;

const double infper = 50.0;
 
int Min(int A, int B) {
	int Result = (A < B) ? A : B;
	if ((A < 0) && (B >= 0)) Result = B;
	if ((B < 0) && (A >= 0)) Result = A;
	if ((A < 0) && (B < 0)) Result = -1;
	return Result;
}

void printM(int *pMatrix, int RowCount, int ColCount) {
	for (int i = 0; i < RowCount; i++) {
		for (int j = 0; j < ColCount; j++)
			printf("%7d", pMatrix[i * ColCount + j]);
		printf("\n");
	}
}

void initializeR(int *pMatrix, int Size) {
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

void term(int *pMatrix) {
	delete[]pMatrix;
}

void initializeP(int *&pMatrix, int& Size) {
	do {
		printf("Enter the number of vertices: ");
		scanf_s("%d", &Size);
		if (Size <= 0)
			printf("The number of vertices should be greater then zero\n");
	} while (Size <= 0);
	pMatrix = new int[Size * Size];
	initializeR(pMatrix, Size);
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

int main(int argc, char* argv[]) {
	int *pMatrix; 
	int Size; 
	time_t start, finish;
	double duration = 0.0;
	printf("SerialFloyd programm\n");
	initializeP(pMatrix, Size);
	printf("Initial matrix:\n");
	printM(pMatrix, Size, Size);
	start = clock();
	SerialFloyd(pMatrix, Size);
	finish = clock();
	printf("\n");
	printf("The matrix after Floyd algorithm:\n");
	printM(pMatrix, Size, Size);
	duration = (finish - start) / double(CLOCKS_PER_SEC);
	printf("Time of execution: %f\n", duration);
	term(pMatrix);
	system("pause");
	return 0;
}
