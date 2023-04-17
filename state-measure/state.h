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
	struct StringHasher {
		std::hash<std::string> hasher;
		size_t operator()(const std::string& t) const {
			return hasher(t);
		}
	};

	template<typename T>
	class StateBase : public std::unordered_map<const std::string, std::shared_ptr<DefaultSchema<T>>, StringHasher> {
		StateBase& operator = (const StateBase&) = delete;
		StateBase& operator = (StateBase&&) noexcept = delete;
		StateBase(const StateBase&) = delete;
		StateBase(StateBase&&) noexcept = delete;
	public:
		StateBase() {};
	};


	template<typename T>
	class State : public std::shared_ptr<StateBase<T>>, public AbstractSchema<T> {
		Basis<T> _basis;
		T _dt;
		template<template<class> class P, class T>
		std::shared_ptr<DefaultSchema<T>>& link(P<T>&& obj) {
			const auto& name = obj->name();
			if ((*this)->find(name) == (*this)->end()) {
				(*this)->insert({ name, std::make_unique<P<T>>(std::forward<P<T>>(obj)) });
			};
			return (*(*this))[name];
		}
	protected:
		StateBase<T>& state;
		template<template<class> class Q, class T>
		std::shared_ptr<DefaultSchema<T>>& add(numerical_schema::type_schema type_schema = numerical_schema::DEFAULT_NUMERICAL_SCHEMA) {
			return this->link(
				numerical_schema::DefaultSchema<T>(
					type_schema,
					std::make_unique<Q<T>>(Q<T>(*this)))
			);
		}

		template<template<class> class Q, class T, size_t N>
		std::shared_ptr<DefaultSchema<T>>& add(numerical_schema::type_schema type_schema = numerical_schema::DEFAULT_NUMERICAL_SCHEMA) {
			this->link(
				numerical_schema::DefaultSchema<T>(
					type_schema,
					std::make_unique<Q<T>>(Q<T>(*this, N)))
			);
		}
	public:
		State(Basis<T>&& basis) : 
			std::shared_ptr<StateBase<T>>(std::make_shared<StateBase<T>>()), 
			_basis(std::move(basis)),
			_dt(0), state(*(*this)) {
		};

		const Basis<T>& basis() {
			return _basis;
		}
		virtual void init() override {
			for (auto& obj : **this) {
				obj.second->init();
			}
		};

		virtual void finalize() override {
			for (auto& obj : **this) {
				obj.second->finalize();
			}
			AbstractSchema<T>::_t += _dt;
		};

		template<typename T>
		friend std::ostream& operator<< (std::ostream& o, const State<T>& b);
	};

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const State<T>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		out << "Measures : \n";
		for (const auto& obj : *b) {
			const auto& value = (*obj.second)->value();
			const auto& rate = (*obj.second)->rate();
			out << " - " << obj.first << ": value = " << value << ", rate = " << rate << std::endl;
		}
		return out;
	};
}