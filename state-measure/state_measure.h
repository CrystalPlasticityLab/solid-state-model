#pragma once
#include "../tensor-matrix/tensor/object.h"
#include <unordered_map>

namespace state {
	template<typename T>
	class State;

	template<typename T>
	class Measure {
		tens::object<T> _rate;
		tens::object<T> _value;
		std::weak_ptr<const State<T>> _state;
		std::string _name;
	public:
		std::string name() { return _name; };
		// -- TODO: add move semantic 
		Measure(const tens::object<T>& obj, std::string name) :
			_name(name),
			_value(obj),
			_rate(obj.get_comp_ref().dim(), obj.get_comp_ref().rank(), tens::FILL_TYPE::ZERO, obj.get_basis_ref()) {
		};

		Measure(tens::object<T>&& obj, std::string name) :
			_name(name),
			_value(std::move(obj)),
			_rate(obj.get_comp_ref().dim(), obj.get_comp_ref().rank(), tens::FILL_TYPE::ZERO, obj.get_basis_ref()) {
		};

		virtual void calc_rate() {};
		virtual void calc_value() {};

		void set_state(const std::shared_ptr<State<T>>& state) {
			_state = state;
		}

		void change_basis(const tens::object<T>& obj) {
			const auto basis = obj.get_basis_ref();
			change_basis(basis);
		}

		void change_basis(const Basis<T>& basis) {
			_rate.change_basis(basis);
			_value.change_basis(basis);
		}

		void add_to_value(const tens::object<T>& obj) {
			_value += obj;
		}

		void add_to_rate(const tens::object<T>& obj) {
			_rate += obj;
		}

		void set_rate(const tens::object<T>& rate) {
			_rate = rate;
		};

		void set_value(const tens::object<T>& value) {
			_value = value;
		};

		tens::object<T> rate() const {
			return _rate;
		};

		tens::object<T> value() const {
			return _value;
		};

		const tens::object<T>& rate_ref() const {
			return _rate;
		};

		const tens::object<T>& value_ref() const {
			return _value;
		};
	};

	struct StringHasher {
		std::hash<std::string> hasher;
		size_t operator()(const std::string& t) const {
			return hasher(t);
		}
	};


	template<typename T>
	class State : public std::unordered_map<const std::string, std::unique_ptr<Measure<T>>, StringHasher> {
		Basis<T> _basis;
		State(Basis<T>&& basis) {
			_basis = std::move(basis);
		}

		template<typename P>
		static void _insert(const std::shared_ptr<State<T>>& state, P&& obj) {
			const auto& name = obj.name();
			if (state->find(name) != state->end()) {
				throw ErrorAccess::Exists();
			};
			obj.change_basis(state->_basis);
			obj.set_state(state);
			state->insert({ name, std::make_unique<Measure<T>>(std::move(obj)) });
		}

		State& operator = (const State&) = delete;
		State& operator = (State&&) noexcept = delete; 
		State(const State&) = delete;
		State(State&&) noexcept = delete;
	public:
		template<typename T>
		friend std::ostream& operator<< (std::ostream& o, const State<T>& b);

		template<typename T>
		friend [[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis);

		template<typename T, typename P>
		friend void insert (const std::shared_ptr<State<T>>& state, P&& obj);
	};

	template<typename T>
	[[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis) {
		const auto obj = new State<T>(std::move(basis));
		return std::shared_ptr<State<T>>(obj);
	}

	template<typename T, typename P>
	void insert(const std::shared_ptr<State<T>>& state, P&& obj){
		state->_insert(state, obj);
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const State<T>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		out << "Measures : \n";
		for (const auto& obj : b){
			const auto& value = obj.second->value_ref().get_comp_ref();
			const auto& rate = obj.second->rate_ref().get_comp_ref();
			out << " - " << obj.first << " (dim=" << value.dim() << ", rank=" << value.rank() << "), value = " << value << ", rate = " << rate << std::endl;
		}
		return out;
	};
}