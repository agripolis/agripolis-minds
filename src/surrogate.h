/*************************************************************************
* This file is part of AgriPoliS-MINDS
*
* AgriPoliS: An Agricultural Policy Simulator
*
* Copyright (c) 2023 Alfons Balmann, Kathrin Happe, Konrad Kellermann et al.
* (cf. AUTHORS.md) at Leibniz Institute of Agricultural Development in
* Transition Economies
*
* SPDX-License-Identifier: MIT
**************************************************************************/


#ifndef SurrogateH
#define SurrogateH
#pragma once

#include <string>
#include <vector>
#include "cppflow/model.h"

using namespace std;
extern shared_ptr<cppflow::model> surrogate_model;
bool dimension_check(string model_dir, int input_dim, int output_dim);
bool predict(vector<float> input, vector<float> &output, string iname, string oname);

#endif