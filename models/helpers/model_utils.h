#pragma once
#include "../../state-measure/state.h"
#include "../relation.h"

namespace model {
	template<class T>
	T parse_json_value(const std::string& name, const json& params) {
		try {
			return params[name].get<T>();
		}
		catch (const std::exception& e) {
			throw std::invalid_argument("Param of plastic model '" + name + "' was not found or has wrong format. Reason " + std::string(e.what()));
		}
	};

	template<typename T>
	class Curve : public std::vector<std::pair<T, T>> {
	private:
		size_t search_lower_bound(T value) const {
			const auto it_begin = this->begin();
			const auto it_end = this->end();
			auto it_last = it_end - 1;
			if (value < it_begin->first || value > it_last->first) {
				throw std::out_of_range("Curve plasticity out of range. Value = " + std::to_string(value));
			}
			auto iter = std::lower_bound(it_begin, it_end, std::make_pair(value, 0.0),
				[](const std::pair<T, T>& a, const std::pair<T, T>& b) {
				return a.first < b.first;
			});
			if (iter == it_begin) iter++;
			if (iter == it_end) iter--;
			return std::distance(it_begin, iter);
		}
	public:
		Curve(const std::vector<std::pair<T, T>>& arr) :
			std::vector<std::pair<T, T>>(arr) {
			if (arr.size() <= 1) {
				throw std::invalid_argument("Curve plasticity: array must contain at least to points");
			}

			for (std::size_t i = 1; i < arr.size(); ++i) {
				if (arr[i].first < arr[i - 1].first) {
					throw std::invalid_argument("Curve plasticity: array must be monotonously increasing by the first component");
				}
			}
		}

		T value(T x) const {
			const size_t pos = search_lower_bound(x); // return pos = [0, last]
			T x_r = (*this)[pos].first;
			T x_l = (*this)[pos - 1].first;
			T y_r = (*this)[pos].second;
			T y_l = (*this)[pos - 1].second;

			return y_l + (y_r - y_l) * (x - x_l) / (x_r - x_l);
		}

		T derivative(T x) const {
			const size_t pos = search_lower_bound(x); // return pos = [0, last]
			T x_r = (*this)[pos].first;
			T x_l = (*this)[pos - 1].first;
			T y_r = (*this)[pos].second;
			T y_l = (*this)[pos - 1].second;

			return (y_r - y_l) / (x_r - x_l);
		}
	};
}