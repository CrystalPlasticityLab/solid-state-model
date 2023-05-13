#pragma once
#include "../state-measure/state.h"
#include "../tensor/object.h"
#include "./relation.h"
#include "./elasticity.h"
#include "./plasticity.h"

namespace model {
	using namespace state;
	using namespace measure;

	template<
		template<template<class> class, template<class> class, class> class Model>
	class ModelFactory {
	public:
		static void test() {};

		template<
			template<class T> class StrainMeasure,
			template<class T> class StressMeasure,
			class T = double>
		static std::shared_ptr<Model<StrainMeasure, StressMeasure, T>> create(const std::string& param_json_file, measure::type_schema type) {
			json params;
			std::ifstream filematerial(param_json_file);
			std::string jsonString;

			try {
				if (filematerial.is_open()) {
					try {
						jsonString = std::string((std::istreambuf_iterator<char>(filematerial)), std::istreambuf_iterator<char>());
						try {
							params = json::parse(jsonString);
						}
						catch (const std::exception& e) {
							throw std::ifstream::failure("Failed during parsing json: " + param_json_file + " \n Reason: " + std::string(e.what()));
						}
					}
					catch (const std::exception& e) {
						filematerial.close();
						throw std::ifstream::failure("Failed during reading file: " + param_json_file + " \n Reason: " + std::string(e.what()));
					}
				}
				else {
					throw std::ios_base::failure("Failed to open file: " + param_json_file);
				}
				try {
					auto result = std::make_shared<Model<StrainMeasure, StressMeasure, T>>(params, measure::type_schema::RATE_CALCULATE);
					return result;
				}
				catch (const std::exception& e) {
					throw std::runtime_error("Model creation was failed. Reason: " + std::string(e.what()));
				}
			} catch (const std::exception& e) {
				std::cout << "Error during model creation: " << e.what();
				exit(1);
			}
		}
	};
};