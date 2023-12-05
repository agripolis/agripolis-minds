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


#include "cppflow/cppflow.h"
#include "cppflow/ops.h"
#include "cppflow/model.h"
#include <iomanip>
#include "surrogate.h"

using namespace std;

bool dimension_check(string model_dir, int input_dim, int output_dim) {
	cout << "VERSION: " << cppflow::version() << endl;

	cppflow::model model(model_dir);
	
	//cout << "size: " << model.get_operations().size() << endl;

	return true;
};


bool predict(vector<float> input, vector<float>& res, string inputname, string outputname) {
	//cout << endl << "VERSION: " << cppflow::version() << endl;
	//_putenv_s("TF_CPP_MIN_LOG_LEVEL", "2");

	//cppflow::model model(model_dir);

	int dim_in = input.size();
	int dim_res = res.size();
	vector<int64_t> ashape = { 1, dim_in };
		
	auto tf_input = cppflow::tensor(input, ashape);
	auto values_in = tf_input.get_data<float>();
	
	auto tf_output = (*surrogate_model)({{inputname, tf_input}}, {outputname});
	//"serving_default_dense_1_input:0", tf_input} }, { "StatefulPartitionedCall:0" });

	auto values = tf_output[0].get_data<float>();
	int cnt = 0;
	for (auto v : values) {
		res[cnt++] = v;
	}
	return true;
};
