#pragma once
#include "container.h"

namespace tens {
	enum class QUATFORM
	{
		REALIMAGE,
		ANGLEAXIS
	};

	template<typename T>
	class quat : public container<T>
	{
	public:
		quat() : container<T>(4, 1, (T)0) { (*this)[0] = (T)1; };
		explicit quat(const quat& q) ;
		explicit quat(const container<T>& q) : container<T>(q) {};
		explicit quat(const T& re, const container<T>& im, QUATFORM type);

		void set_im(const container<T>& im);
		inline container<T> get_im() const;
		inline T re() const;

		inline quat& operator = (const container<T>& v);
		inline quat operator * (const quat& rhs) const;
		inline quat operator * (const T& mult) const;
		inline quat operator * () const;
		inline quat operator ! () const;

		static friend container<T> get_ort_matrix(const quat<T>& q);
		//static inline matrix<T, 3> get_ort_matrix(const T& angle, const container<T, 3, 1>& axis);
	private:
	};

	template<typename T>
	quat<T>::quat(const T& re, const container<T>& im, QUATFORM type) {
		T angleto2 = (T)0;
		switch (type)
		{
		case QUATFORM::REALIMAGE:
			this->set(0, re);
			set_im(im);
			break;
		case QUATFORM::ANGLEAXIS:
			angleto2 = re * (T)0.5;
			this->set(0, cos(angleto2));
			set_im(im.get_normalize() * sin(angleto2));
			break;
		default:
			break;
		}
	};

	template<typename T>
	quat<T>::quat(const quat& q ) {
		this->set(0, q.get(0));
		set_im(q.get_im());
	};

	template<typename T>
	inline void quat<T>::set_im(const container<T>& im){
		(*this)[1] = im[0];
		(*this)[2] = im[1];
		(*this)[3] = im[2];
	}

	template<typename T>
	inline container<T> quat<T>::get_im() const{
		container<T, 3, 1> res;
		const std::array<T, 4>& v = (*this)();
		res[0] = v[1];
		res[1] = v[2];
		res[2] = v[3];
		return res;
	}

	template<typename T>
	inline T quat<T>::re() const{
		return this->get(0);
	}

	template<typename T>
	inline quat<T>& quat<T>::operator = (const container<T>& v){
		static_cast<quat<T>&>(container<T, 4, 1>::operator=(v));
		return *this;
	}

	template<typename T>
	inline quat<T> quat<T>::operator * () const{
		quat<T> res(*this);
		res.set_im(-this->get_im());
		return res;
	}

	template<typename T>
	inline quat<T> quat<T>::operator ! () const{
		quat<T> res(*(*this));
		T norm2 = res * res;
		return res / norm2;
	}

	template<typename T>inline quat<T> quat<T>::operator * (const T& mult) const{
		return static_cast<quat<T>>(container<T, 4, 1>::operator*(mult));
	}

	template<typename T>
	inline quat<T>  quat<T>::operator * (const quat<T>& rhs) const{
		quat<T> res;
		container<T, 3, 1> lv = get_im();
		container<T, 3, 1> rv = rhs.get_im();
		T lr = re();
		T rr = rhs.re();

		res.set(0, lr * rr - lv * rv);
		res.set_im(rv * lr + lv * rr + lv.vector_product(rv));
		return res;
	}

	template<typename T>
	static container<T> get_ort_matrix(const quat<T>& _q){
		T wx, wy, wz, xx, yy, yz, xy, xz, zz;
		std::array<std::array<T, 3>, 3> m;
		container<T> q = get_normalize(_q);

		xx = 2.0 * q[1] * q[1];   xy = 2.0 * q[1] * q[2];   xz = 2.0 * q[1] * q[3];
		yy = 2.0 * q[2] * q[2];   yz = 2.0 * q[2] * q[3];   zz = 2.0 * q[3] * q[3];
		wx = 2.0 * q[0] * q[1];   wy = 2.0 * q[0] * q[2];   wz = 2.0 * q[0] * q[3];

		m[0][0] = 1.0 - (yy + zz); /*+*/ m[0][1] = xy - wz;               /*+*/  m[0][2] = xz + wy;				 /*+*/
		m[1][0] = xy + wz;         /*+*/ m[1][1] = 1.0 - (xx + zz);       /*+*/	 m[1][2] = yz - wx;				 /*+*/
		m[2][0] = xz - wy;         /*+*/ m[2][1] = yz + wx;		          /*+*/  m[2][2] = 1.0 - (xx + yy);/*+*/
		return container<T>(Matrix<double, 3>(m));
	}

	//template<typename T>
	//inline matrix<T, 3> quat<T>::get_ort_matrix(const T& angle, const container<T, 3, 1>& axis){
	//	return quat<T>(angle, axis, QUATFORM::ANGLEAXIS).get_ort_matrix();
	//}
};