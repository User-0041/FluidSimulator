#include <iostream>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "Simulation.h"
void drawPressure(cv::Mat& image, std::vector<float>* p, int numX, int numY) {
    for (int i = 1; i < numX; i++) {
        for (int j = 1; j < numY; j++) {
            int x = i * 5; // Scale factor for visualization
            int y = j * 5;

            // Map pressure to a color gradient (blue to red)
            int blue = static_cast<int>(255 * (1 - p->at(i * numY + j) / 100.0)); // Adjust the scale for pressure
            int red = static_cast<int>(255 * (p->at(i * numY + j) / 100.0));

            cv::Scalar pressureColor(blue, 0, red);
            cv::rectangle(image, cv::Rect(x, y, 5, 5), pressureColor, cv::FILLED);
        }
    }
}

void drawDensity(cv::Mat& image, std::vector<float>* m, int numX, int numY) {
    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            int x = i * 5; // Scale factor for visualization
            int y = j * 5;

            // Map density to a color gradient (blue to red)
            int blue = static_cast<int>(255 * (1 - m->at(i * numY + j) / 100.0)); // Adjust the scale for density
            int red = static_cast<int>(255 * (m->at(i * numY + j) / 100.0));

            cv::Scalar densityColor(blue, 0, red);
            cv::rectangle(image, cv::Rect(x, y, 5, 5), densityColor, cv::FILLED);
        }
    }
}

void drawArrows(cv::Mat& image, std::vector<float>* u, std::vector<float>* v, int numX, int numY) {
    float arrowScale = 10.0; // Adjust the scale for arrow length

    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            int x = i * 5; // Scale factor for visualization
            int y = j * 5;

            int arrowX = static_cast<int>(x + arrowScale * u->at(i * numY + j));
            int arrowY = static_cast<int>(y + arrowScale * v->at(i * numY + j));

            // Calculate velocity magnitude
            float velocityMagnitude = sqrt(u->at(i * numY + j) * u->at(i * numY + j) + v->at(i * numY + j) * v->at(i * numY + j));

            // Map velocity magnitude to a color gradient (blue to red)
            int blue = static_cast<int>(255 * (1 - velocityMagnitude / arrowScale));
            int red = static_cast<int>(255 * (velocityMagnitude / arrowScale));

            cv::Scalar arrowColor(blue, 0, red);
            cv::arrowedLine(image, cv::Point(x, y), cv::Point(arrowX, arrowY), arrowColor, 1, 8, 0, 0.1);
        }
    }
}
void main() {
	Simulation* s = new Simulation(1000.0f, 200, 200, 1);
    s->initializeS();
	float dt = 1.0f / 30.0f;
    float g = -1.0;
    // Create an OpenCV window
    cv::namedWindow("Fluid Simulation", cv::WINDOW_NORMAL);
    cv::resizeWindow("Fluid Simulation", s->numX * 5, s->numY * 5);
    cv::namedWindow("Fluid Pressure", cv::WINDOW_NORMAL);
    cv::resizeWindow("Fluid Pressure", s->numX * 5, s->numY * 5);
    cv::namedWindow("Density Visualization", cv::WINDOW_NORMAL);
    cv::resizeWindow("Density Visualization", s->numX * 5, s->numY * 5);



	while (true)
	{
		s->integrate(dt,g);
		s->solveIncompressibility(100, dt);
        s->extrapolate();
		s->advectVel(dt);
		s->advectSmoke(dt);


        cv::Mat velocityVisualization(s->numX * 5, s->numY * 5, CV_8UC3, cv::Scalar(255, 255, 255));
        cv::Mat pressureVisualization(s->numX * 5, s->numY * 5, CV_8UC3, cv::Scalar(255, 255, 255));
        cv::Mat densityVisualization(s->numX * 5, s->numY * 5, CV_8UC3, cv::Scalar(255, 255, 255));

        drawArrows(velocityVisualization, s->u, s->v, s->numX , s->numY );
        drawPressure(pressureVisualization, s->p, s->numX , s->numY );
        drawDensity(densityVisualization, s->m, s->numX , s->numY );


        cv::imshow("Fluid Simulation", velocityVisualization);
        cv::imshow("Fluid Pressure", pressureVisualization);
        cv::imshow("Density Visualization", densityVisualization);
        
		char key = cv::waitKey(1);
		if (key == 27) {  // ASCII value for the Escape key
			break;
		}

	}
}