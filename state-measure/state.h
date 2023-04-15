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
	class StateBase : public std::unordered_map<const std::string, std::unique_ptr<DefaultSchema<T>>, StringHasher>, public AbstractSchema<T> {
		Basis<T> _basis;
		T _dt;
		T _t;

		template<template<class> class P, class T>
		void link(P<T>&& obj) {
			const auto& name = obj->name();
			if (this->find(name) != this->end()) {
				throw ErrorAccess::Exists();
			};
			this->insert({ name, std::make_unique<P<T>>(std::forward<P<T>>(obj)) });
		}

		StateBase& operator = (const StateBase&) = delete;
		StateBase& operator = (StateBase&&) noexcept = delete;
		StateBase(const StateBase&) = delete;
		StateBase(StateBase&&) noexcept = delete;
	public:
		StateBase(Basis<T>&& basis) noexcept {
			_dt = T();
			_t = T();
			_basis = std::move(basis);
		}
		const Basis<T>& basis() {
			return _basis;
		}
		T t() { return _t; };
		template<typename T>
		friend std::ostream& operator<< (std::ostream& o, const StateBase<T>& b);

		virtual void init() override {
			for (auto& obj : *this) {
				obj.second->init();
			}
		};

		virtual void calc(T dt) override {
			for (auto& obj : *this) {
				obj.second->calc(dt);
			}
			_dt = dt;
		};

		virtual void finalize() override {
			for (auto& obj : *this) {
				obj.second->finalize();
			}
			_t += _dt;
		};
	};


	template<typename T>
	std::ostream& operator<<(std::ostream& out, const StateBase<T>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		out << "Measures : \n";
		for (const auto& obj : b) {
			const auto& value = (*obj.second)->value();
			const auto& rate = (*obj.second)->rate();
			out << " - " << obj.first << ": value = " << value << ", rate = " << rate << std::endl;
		}
		return out;
	};

	template<typename T>
	class State : public std::shared_ptr<StateBase<T>> {
		template<template<class> class P, class T>
		void link(P<T>&& obj) {
			const auto& name = obj->name();
			if ((*this)->find(name) != (*this)->end()) {
				throw ErrorAccess::Exists();
			};
			(*this)->insert({ name, std::make_unique<P<T>>(std::forward<P<T>>(obj)) });
		}
	protected:
		template<template<class> class Q, class T>
		void add(numerical_schema::type_schema type_schema) {
			this->link(
				numerical_schema::DefaultSchema<T>(
					type_schema,
					std::make_unique<Q<T>>(Q<T>(*this)))
			);
		}

		template<template<class> class Q, class T, size_t N>
		void add(numerical_schema::type_schema type_schema) {
			this->link(
				numerical_schema::DefaultSchema<T>(
					type_schema,
					std::make_unique<Q<T>>(Q<T>(*this, N)))
			);
		}
	public:
		State(Basis<T>&& basis) : std::shared_ptr<StateBase<T>>(std::make_shared<StateBase<T>>(std::move(basis))) {};
		~State() {
			return;
		}
	};
}