#pragma once
#include "../tensor-matrix/tensor/object.h"
#include "../tensor-matrix/state-measure/measure.h"
#include "../tensor-matrix/state-measure/schema.h"
#include "../tensor-matrix/state-measure/strain.h"
#include "../tensor-matrix/state-measure/stress.h"
#include "../tensor-matrix/state-measure/scalar.h"
#include <unordered_map>

namespace state {
	using namespace measure;
	using namespace numerical_schema;

	/*
		Base class of Material Point with basis
	*/
	template<class T, size_t DIM>
	class MaterialPoint : public AbstractSchema<T> {
		Basis<T, DIM> _basis;
	protected:
		std::shared_ptr<Json::Value> _params;
		virtual void parse_json_params(const Json::Value& params) {};
	public:
		virtual void init() override {};
		virtual void calc(T dt) override {};
		virtual void finalize()  override {};
		MaterialPoint(const Json::Value& params, measure::type_schema type) :
			_basis(tens::create_basis<T, 3>(tens::DEFAULT_ORTH_BASIS::RANDOM))
		{
		};
		const Basis<T, DIM>& basis() {
			return _basis;
		}

		const std::shared_ptr<const Json::Value>& param() const {
			return _params;
		}
		template<typename T, size_t DIM>
		friend std::ostream& operator<< (std::ostream& o, const MaterialPoint<T, DIM>& b);

		virtual std::ostream& print_measures(std::ostream& o) const = 0;
	};

	template<typename T, size_t DIM>
	std::ostream& operator<<(std::ostream& out, const MaterialPoint<T, DIM>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		b.print_measures(out);
		return out;
	};
}