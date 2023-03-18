#pragma once
#include "../tensor-matrix/tensor/object.h"
#include "../tensor-matrix/state-measure/measure.h"
#include <unordered_map>

namespace state {
	using namespace measure;
	struct StringHasher {
		std::hash<std::string> hasher;
		size_t operator()(const std::string& t) const {
			return hasher(t);
		}
	};

	template<typename T>
	class State : public std::unordered_map<const std::string, std::unique_ptr<Measure<T>>, StringHasher> {
		Basis<T> _basis;
		T _dt;

		State(Basis<T>&& basis) {
			_dt = 0.123456;
			_basis = std::move(basis);
		}

		template<typename P>
		void _insert(P&& obj) {
			const auto& name = obj.name();
			if (this->find(name) != this->end()) {
				throw ErrorAccess::Exists();
			};
			this->insert({ name, std::make_unique<Measure<T>>(std::forward<P>(obj)) });
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
		template<typename T>
		friend std::ostream& operator<< (std::ostream& o, const State<T>& b);

		template<typename T>
		friend [[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis);

		template<typename T, typename P>
		friend void insert (const std::shared_ptr<State<T>>& state, P&& obj);
	};

	template<typename T>
	[[nodiscard]] std::shared_ptr<State<T>> create(Basis<T>&& basis) {
		return std::shared_ptr<State<T>>(new State<T>(std::move(basis)));
	}

	template<typename T, typename P>
	void insert(const std::shared_ptr<State<T>>& state, P&& obj){
		state->_insert(obj);
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const State<T>& b) {
		out << "Basis    : " << *b._basis << std::endl;
		out << "Measures : \n";
		for (const auto& obj : b){
			const auto& value = obj.second->value().get_comp_ref();
			const auto& rate = obj.second->rate().get_comp_ref();
			out << " - " << obj.first << " (dim=" << value.dim() << ", rank=" << value.rank() << "), value = " << value << ", rate = " << rate << std::endl;
		}
		return out;
	};
}