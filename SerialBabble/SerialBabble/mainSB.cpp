#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
 
using namespace std;

const double RandomMultiplier = 1000.0;

void RandomData(double *&pData, int& DataSize) {
	srand((unsigned)time(0));
	for (int i = 0; i < DataSize; i++)
		pData[i] = double(rand()) / RAND_MAX * RandomMultiplier;
}

void initialize(double *&pData, int& DataSize) {
	do {
		printf("Enter the size of data: ");
		scanf_s("%d", &DataSize);
		if (DataSize <= 0)
			printf("Data size should be greater than zero\n");
	} while (DataSize <= 0);
	printf("Sorting %d data items\n", DataSize);
	pData = new double[DataSize];
	RandomData(pData, DataSize);
}

void term(double *pData) {
	delete[]pData;
}

void SerialBubble(double *pData, int DataSize) {
	double Tmp;
	for (int i = 1; i < DataSize; i++)
		for (int j = 0; j < DataSize - i; j++)
			if (pData[j] > pData[j + 1]) {
				Tmp = pData[j];
				pData[j] = pData[j + 1];
				pData[j + 1] = Tmp;
			}
}

void PrintData(double *pData, int DataSize) {
	for (int i = 0; i < DataSize; i++)
		printf("%7.4f ", pData[i]);
	printf("\n");
}

int main(int argc, char *argv[]) {
	double *pData = 0;
	int DataSize = 0;
	time_t start, finish;
	double duration = 0.0;

	printf("Serial bubble sort program\n");
	 
	initialize(pData, DataSize);
	cout << "************************" << endl;
	printf("Data before sorting:\n");
	PrintData(pData, DataSize);
	start = clock();
	SerialBubble(pData, DataSize);
	finish = clock();
	cout << "***********************" << endl;
	printf("Data after sorting\n");
	PrintData(pData, DataSize);
	duration = (finish - start) / double(CLOCKS_PER_SEC);
	cout << "***********************" << endl;
	printf("Time of execution: %f\n", duration);
	term(pData);

	system("pause");

	return 0;
}