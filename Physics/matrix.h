#pragma once
#include <vector>
#include <math.h>
#include <bits/stdc++.h>

float getDeterminant(const std::vector<std::vector<float>> vect);
std::vector<std::vector<float>> getTranspose(const std::vector<std::vector<float>> matrix1);
std::vector<std::vector<float>> getCofactor(const std::vector<std::vector<float>> vect);
std::vector<std::vector<float>> getInverse(const std::vector<std::vector<float>> vect);

std::vector<std::vector<float>> zeros(int width, int height);
float det(std::vector<std::vector<float>> A);
float determinantOfMatrix(std::vector<std::vector<float>> mat, int n);
std::vector<float> linspace(int startX, int endX, int steps);