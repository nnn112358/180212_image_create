#include <iostream> 
#include <cmath>
#include <string>
#include <sstream>

#include <vector>
#include <algorithm>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <ctype.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef _DEBUG
#define CV_EXT "d.lib"
#else
#define CV_EXT ".lib"
#endif
#pragma comment(lib, "opencv_world320" CV_EXT) // OpenCV3.3.0�̏ꍇ�́A"opencv_core330"��

using namespace std;
using namespace cv;

const double img_dpi = 9600;
const double um2dot = img_dpi / 25.4 / 1000.0;
const double mm2dot = img_dpi / 25.4 ;
const double dot2um = 1.0/um2dot;
const double dot2mm = 1.0 / mm2dot;


void img_info(Mat &img);
int int_random(int max);
double uniform_random(double min, double max);
double gaussian_random(double mean, double std);


int main(int argc, char* argv[]){

	double img_width_mm = 32.0;
	double img_height_mm = 130;

	double gray_white = 20;			//White
	double gray_black_ave = 228;	//Black
	double gray_black_stdev = 0;

	double dot_r_ave = 80 * um2dot;
	double dot_r_stdev = 0.0*um2dot;

	double yy_sin_amp = 3.0*um2dot;
	double yy_sin_pitch = 1.0;
	double xx_sin_amp = 0.0*um2dot;
	double xx_sin_pitch = 10.0;

	//�摜�̍쐬
	int img_width = img_width_mm* mm2dot;
	int img_height = img_height_mm * mm2dot;

	Mat img_proc= Mat::zeros( img_height, img_width, CV_8U);
	
	Mat img_proc_200dpi = Mat::zeros(img_height*200.0 / img_dpi, img_width*200.0/ img_dpi, CV_8U);
	Mat img_proc_400dpi = Mat::zeros(img_height*400.0 / img_dpi, img_width*400.0/ img_dpi, CV_8U);
	Mat img_proc_1200dpi = Mat::zeros(img_height*1200.0 / img_dpi, img_width*1200.0 / img_dpi, CV_8U);
	Mat img_proc_2400dpi = Mat::zeros(img_height*2400.0 / img_dpi, img_width*2400.0 / img_dpi, CV_8U);

	img_proc += 255-gray_white;
	img_info(img_proc);

	std::vector<double> screen_x_um;
	std::vector<double> screen_y_um;
	std::vector<double> screen_r_um;
	std::vector<double> screen_r_dense;

	double screen_grid_um= 25.4  * 1000.0 * 16.0 / 2400.0;
	cout << "screen_grid_num" << endl;


	for (long i = 0; i < img_width / 16; i++) {
		for (long j = 0; j < img_height/16; j++) {
			double xx = screen_grid_um*i;
			double yy = screen_grid_um*j;
			screen_x_um.push_back(xx);
			screen_y_um.push_back(yy);
			screen_x_um.push_back(xx + screen_grid_um / 2.0);
			screen_y_um.push_back(yy + screen_grid_um / 2.0);

			double dot_r1 = gaussian_random(dot_r_ave/2.0, dot_r_stdev);
			double dot_r2 = gaussian_random(dot_r_ave / 2.0, dot_r_stdev);
			if (dot_r1 < 0) dot_r1 = 0.0;
			if (dot_r2 < 0) dot_r2 = 0.0;
			
			screen_r_um.push_back(dot_r1);
			screen_r_um.push_back(dot_r2);

			double dot_dense1 = gaussian_random(gray_black_ave, gray_black_stdev);
			double dot_dense2 = gaussian_random(gray_black_ave, gray_black_stdev);
			if (dot_dense1 < 0)dot_dense1 = 0;
			if (dot_dense1 > 255)dot_dense1 = 255;
			if (dot_dense2 < 0)dot_dense2 = 0;
			if (dot_dense2 > 255)dot_dense2 = 255;

			screen_r_dense.push_back(dot_dense1);
			screen_r_dense.push_back(dot_dense2);
		}
	}

	cout << "circle_generation" << endl;


	for (long i = 0; i < screen_x_um.size(); i++) {
			double xx_um = screen_x_um[i];
			double yy_um = screen_y_um[i];

			double xx_delt = xx_sin_amp*sin(2.0*M_PI*yy_um / 1000.0 / xx_sin_pitch);
			xx_um += xx_delt;

			double yy_delt = yy_sin_amp*sin(2.0*M_PI*yy_um / 1000.0 / yy_sin_pitch);
			yy_um += yy_delt;

			double dot_r = screen_r_um[i];
			int dot_dense =255- screen_r_dense[i];
			cv::circle(img_proc, cv::Point(xx_um*um2dot, yy_um*um2dot), dot_r, cv::Scalar(dot_dense), -1, CV_AA);

			if (i == screen_x_um.size()*0.25) {
				cout << "Progress 25%" << endl;
			}
			else if (i == screen_x_um.size()*0.5) {
				cout << "Progress 50%" << endl;
			}
			else if (i == screen_x_um.size()*0.75) {
				cout << "Progress 75%" << endl;
			}
	}
	cout << "resize img_proc_2400dpi" << endl;
	cv::resize(img_proc, img_proc_2400dpi, img_proc_2400dpi.size(), cv::INTER_AREA);

	cout << "resize img_proc_1200dpi" << endl;
	cv::resize(img_proc_2400dpi, img_proc_1200dpi, img_proc_1200dpi.size(), cv::INTER_AREA);

	cout << "resize img_proc_400dpi" << endl;
	cv::resize(img_proc_1200dpi, img_proc_400dpi, img_proc_400dpi.size(), cv::INTER_AREA);

	cout << "resize img_proc_200dpi" << endl;
	cv::resize(img_proc_400dpi, img_proc_200dpi, img_proc_200dpi.size(), cv::INTER_AREA);

	cout << "imwrite" << endl;

	cv::imwrite("output.tif", img_proc);
	cv::imwrite("output_200dpi.tif", img_proc_200dpi);
	cv::imwrite("output_400dpi.tif", img_proc_400dpi);
	cv::imwrite("output_1200dpi.tif", img_proc_1200dpi);

	return 0;
}



void img_info(Mat &img) {

	// �s��
	std::cout << "rows: " << img.rows << std::endl;
	// ��
	std::cout << "cols: " << img.cols << std::endl;


	// �������i�摜�Ȃ̂ŏc�E����2�����j
	std::cout << "dims: " << img.dims << std::endl;
	// �T�C�Y�i2�����̏ꍇ�j
	std::cout << "size[]: " << img.size().width << "," << img.size().height << std::endl;
	// �r�b�g�[�xID
	std::cout << "depth (ID): " << img.depth() << "(=" << CV_8U << ")" << std::endl;
	// �`�����l����
	std::cout << "channels: " << img.channels() << std::endl;
	// �i�����`�����l�����琬��j1�v�f�̃T�C�Y [�o�C�g�P��]
	std::cout << "elemSize: " << img.elemSize() << "[byte]" << std::endl;
	// 1�v�f����1�`�����l�����̃T�C�Y [�o�C�g�P��]
	std::cout << "elemSize1 (elemSize/channels): " << img.elemSize1() << "[byte]" << std::endl;
	// �^�C�v
	std::cout << "type (ID): " << img.type() << "(=" << CV_8UC3 << ")" << std::endl;
	// �v�f�̑���
	std::cout << "total: " << img.total() << std::endl;
	// �X�e�b�v�� [�o�C�g�P��]
	std::cout << "step: " << img.step << "[byte]" << std::endl;
	// 1�X�e�b�v���̃`�����l������
	std::cout << "step1 (step/elemSize1): " << img.step1() << std::endl;
	// �f�[�^�͘A�����H
	std::cout << "isContinuous: " << (img.isContinuous() ? "true" : "false") << std::endl;
	// �����s�񂩁H
	std::cout << "isSubmatrix: " << (img.isSubmatrix() ? "true" : "false") << std::endl;
	// �f�[�^�͋󂩁H
	std::cout << "empty: " << (img.empty() ? "true" : "false") << std::endl;

}



// max���ő�Ƃ���int�^�̗����𐶐�
int int_random(int max) {
	int i;
	return (int)(max*(rand() / (RAND_MAX + 1.0)));
}

// min���ŏ��Cmax���ő�Ƃ���double�^�̗����𐶐�
double uniform_random(double min, double max) {
	return min + (rand() / (double)RAND_MAX) * (max - min);
}

// ����mean�C�W���΍�std�̃K�E�X���z�̗����𐶐�
double gaussian_random(double mean, double std) {
	const double norm = 1.0 / (RAND_MAX + 1.0);
	double u = 1.0 - rand() * norm;
	double v = rand() * norm;
	double z = sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
	return mean + std * z;
}

