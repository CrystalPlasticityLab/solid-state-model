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
		T _t;

		State(Basis<T>&& basis) noexcept {
			_dt = T();
			_basis = std::move(basis);
		}

		template<template<class> class P, class T>
		void link(P<T>&& obj) {
			const auto& name = obj->name();
			if (this->find(name) != this->end()) {
				throw ErrorAccess::Exists();
			};
			this->insert({ name, std::make_unique<P<T>>(std::forward<P<T>>(obj)) });
		}

		State& operator = (const State&) = delete;
		State& operator = (State&&) noexcept = delete;
		State(const State&) = delete;
		State(State&&) noexcept = delete;
	public:
		const Basis<T>& basis() {
			return _basis;
		}
		T dt() { return _dt; };
		T t() { return _t; };
		void set_time_integration_step(T dt) { _dt = dt; };
		template<typename T>
		friend std::ostream& operator<< (std::ostream& o, const State<T>& b);

		template<typename T>
		friend [[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis) noexcept;

		template<template<class> class Q, class T>
		friend void add(std::shared_ptr<State<T>>& state, numerical_schema::type_schema type_schema);

		template<template<class> class Q, class T, size_t N>
		friend void add(std::shared_ptr<State<T>>& state, numerical_schema::type_schema type_schema);

		virtual void init() override {
			for (auto& obj : *this) {
				obj.second->init();
			}
		};

		virtual void calc() override {
			for (auto& obj : *this) {
				obj.second->calc();
			}
		};

		virtual void finalize() override {
			for (auto& obj : *this) {
				obj.second->finalize();
			}
			_t += _dt;
		};
	};

	template<template<class> class Q, class T>
	void add(std::shared_ptr<State<T>>& state, numerical_schema::type_schema type_schema) {
		state->link(
			numerical_schema::DefaultSchema<T>(
				type_schema, 
				std::make_unique<Q<T>>(Q<T>(state)), 
				state->_dt)
		);
	}

	template<template<class> class Q, class T, size_t N>
	void add(std::shared_ptr<State<T>>& state, numerical_schema::type_schema type_schema) {
		state->link(
			numerical_schema::DefaultSchema<T>(
				type_schema,
				std::make_unique<Q<T>>(Q<T>(state, N)),
				state->_dt)
		);
	}

	template<typename T>
	[[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis) noexcept {
		return std::shared_ptr<State<T>>(new State<T>(std::move(basis)));
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const State<T>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		out << "Measures : \n";
		for (const auto& obj : b){
			const auto& value = (*obj.second)->value();
			const auto& rate = (*obj.second)->rate();
			out << " - " << obj.first << ": value = " << value << ", rate = " << rate << std::endl;
		}
		return out;
	};
}