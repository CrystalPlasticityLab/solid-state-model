#pragma once
#include "../tensor-matrix/tensor/object.h"
#include "../tensor-matrix/state-measure/measure.h"
#include "../tensor-matrix/state-measure/schema.h"
#include "../tensor-matrix/state-measure/strain.h"
#include "../tensor-matrix/state-measure/stress.h"
#include <unordered_map>

namespace state {
	using namespace measure;
	using namespace numerical_schema;
	struct StringHasher {
		std::hash<std::string> hasher;
		size_t operator()(const std::string& t) const {
			return hasher(t);
		}
	};

	template<typename T>
	class State : public std::unordered_map<const std::string, std::unique_ptr<DefaultSchema<T>>, StringHasher>, public AbstractSchema {
		Basis<T> _basis;
		T _dt;

		State(Basis<T>&& basis) {
			_dt = 0.123456;
			_basis = std::move(basis);
		}

		State& operator = (const State&) = delete;
		State& operator = (State&&) noexcept = delete;
		State(const State&) = delete;
		State(State&&) noexcept = delete;
	public:
		const Basis<T>& basis() {
			return _basis;
		}
		T& dt() { return _dt; };
		void set_time_integration_step(T dt) { _dt = dt; };
		template<typename T>
		friend std::ostream& operator<< (std::ostream& o, const State<T>& b);

		template<typename T>
		friend [[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis);

		template<typename T>
		friend void link(const std::shared_ptr<State<T>>& state, DefaultSchema<T>&& obj);

		virtual void init() override {
			for (auto& obj : *this) {
				obj.second->init();
			}
		};
		virtual void step() override {
			for (auto& obj : *this) {
				obj.second->step();
			}
		};
		virtual void finalize() override {
			for (auto& obj : *this) {
				obj.second->finalize();
			}
		};

		template<typename T>
		static void add(std::string measure_name, std::shared_ptr<State<T>>& state, numerical_schema::type_schema type_schema);
	};

	template<typename T>
	static void add(std::string measure_name, std::shared_ptr<State<T>>& state, numerical_schema::type_schema type_schema) {
		if (measure_name == strain::DEFORM_GRADIENT) {
			numerical_schema::DefaultSchema<T>(type_schema, strain::GradDeform<T>(state), state);
		}
		else if (measure_name == stress::CAUCHY) {
			numerical_schema::DefaultSchema<T>(type_schema, stress::CaushyStress<T>(state), state);
		}
	}

	template<typename T>
	[[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis) {
		return std::shared_ptr<State<T>>(new State<T>(std::move(basis)));
	}

	template<typename T>
	void link(const std::shared_ptr<State<T>>& state, DefaultSchema<T>&& obj) {
		const auto& name = obj.name();
		if (state->find(name) != state->end()) {
			throw ErrorAccess::Exists();
		};
		state->insert({ name, std::make_unique<DefaultSchema<T>>(std::move(obj))});
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const State<T>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		out << "Measures : \n";
		for (const auto& obj : b){
			const auto& value = obj.second->get_value();
			const auto& rate = obj.second->get_rate();
			out << " - " << obj.first << " (dim=" << value.dim() << ", rank=" << value.rank() << "), value = " << value << ", rate = " << rate << std::endl;
		}
		return out;
	};
}