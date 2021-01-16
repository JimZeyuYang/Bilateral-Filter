#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>

int addrY(int posnH, int posnW);
int addrU(int posnH, int posnW);
int addrV(int posnH, int posnW);
float distance(int x, int y, int i, int j);
double gaussian(float x, double sigma);
void bilateralFilter(int x, int y, int J, int K, unsigned char inputBuffer[], unsigned char outputBuffer[]);
double calculatePSNR(unsigned char inputBuffer[], unsigned char outputBuffer[]);

int width = 832;
int height = 480;
int frames = 10;
int blockSize = 16;

int kernelSize = 5;
int sigmaI = 8;
int sigmaS = 8;

int frameSize = width * height * 1.5;

int main() {

	std::ifstream file_org("BQMall_832x480_60.YUV", std::ios::binary);
	std::ofstream file_copy("NEW.YUV", std::ios::binary);
	
	if (!file_org.is_open()) {
		std::cout << "The input file couldn't be opened" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (!file_copy.is_open()) {
		std::cout << "The output file couldn't be created" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	unsigned char* inputBuffer = new unsigned char[frameSize];
	unsigned char* outputBuffer = new unsigned char[frameSize];

	double* PSNR = new double[frames];

	for (double i = 0; i < frames; i++) {
		file_org.seekg(i*frameSize);
		file_org.read(reinterpret_cast<char*>(inputBuffer), frameSize);
		for (int j = 0; j < height / blockSize; j++) {
			for (int k = 0; k < width / blockSize; k++) {
				int blockH = blockSize;
				int blockW = blockSize;
				if (j == height / blockSize) {
					blockH = height % blockSize;
				}
				if (k == width / blockSize) {
					blockW = width % blockSize;
				}
				for (int l = 0; l < blockH; l++) {
					for (int m = 0; m < blockW; m++) {
						bilateralFilter(l, m, j, k, inputBuffer, outputBuffer);

						outputBuffer[addrU(j * blockSize + l, k * blockSize + m)] 
							= inputBuffer[addrU(j * blockSize + l, k * blockSize + m)];

						outputBuffer[addrV(j * blockSize + l, k * blockSize + m)] 
							= inputBuffer[addrV(j * blockSize + l, k * blockSize + m)];
					}
				}
			}
		}
		PSNR[static_cast<int>(i)] = calculatePSNR(inputBuffer, outputBuffer);

		std::cout << PSNR[static_cast<int>(i)] << " ";

		file_copy.seekp(i * frameSize);
		file_copy.write(reinterpret_cast<char*>(outputBuffer), frameSize);
	}


	delete[] inputBuffer;
	delete[] outputBuffer;

	file_org.close();
	file_copy.close();
	
	return 0;
}

int addrY(int posnH, int posnW) {
	return posnH * width + posnW;
}

int addrU(int posnH, int posnW) {
	return height * width + posnH * width / 4 + posnW / 2;
}

int addrV(int posnH, int posnW) {
	return height * width * 1.25 + posnH * width / 4 + posnW / 2;
}

float distance(int x, int y, int i, int j) {
	return float(sqrt(pow(x - i, 2) + pow(y - j, 2)));
}

double gaussian(float x, double sigma) {
	return exp(-(pow(x, 2)) / (2 * pow(sigma, 2))) / (2 * 3.14159265358 * pow(sigma, 2));
}

void bilateralFilter(int x, int y, int J, int K,unsigned char inputBuffer[], unsigned char outputBuffer[]) {
	double iFiltered = 0;
	double wP = 0;
	int neighborX = 0;
	int neighborY = 0;

	for (int i = 0; i < kernelSize; i++) {
		for (int j = 0; j < kernelSize; j++) {
			neighborX = x - (kernelSize / 2 - i);
			neighborY = y - (kernelSize / 2 - j);

			if (neighborX >= 0 && neighborY >= 0) {
				double gi = gaussian(inputBuffer[addrY(J * blockSize + x, K * blockSize + y)]
					- inputBuffer[addrY(J * blockSize + neighborX, K * blockSize + neighborY)], sigmaI);
				double gs = gaussian(distance(x, y, neighborX, neighborY), sigmaS);
				double w = gi * gs;
				iFiltered += w * inputBuffer[addrY(J * blockSize + neighborX, K * blockSize + neighborY)];
				wP += w; 
			}
		}
	}
	iFiltered /= wP;
	outputBuffer[addrY(J * blockSize + x, K * blockSize + y)] = iFiltered;
}

double calculatePSNR(unsigned char inputBuffer[], unsigned char outputBuffer[]) {
	double MSE = 0;
	for (int i = 0; i <= height; i++) {
		for (int j = 0; j <= width; j++) {
			MSE += (inputBuffer[addrY(i, j)] - outputBuffer[addrY(i, j)]);
		}
	}
	return 10 * log10(pow(155, 2) / (MSE / (static_cast<double>(width) * height)));
}