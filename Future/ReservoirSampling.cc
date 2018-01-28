// Sketches Library
//
// Copyright (C) 2005 Marios Hadjieleftheriou
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Email:
//    mhadji@gmail.com

#include <Sketches.h>

template<class T> Sketches::ReservoirSampling_impl<T>::ReservoirSampling_impl(
	double epsilon,
	double delta,
	Tools::Random& r,
	unsigned long threshold
) : m_sampleSize(static_cast<unsigned long>((0.5 * std::log(2.0 / delta)) / std::pow(epsilon, 2.0))),
    m_t(0),
    m_T(threshold),
    m_pRandom(&r)
{
}

template<class T> Sketches::ReservoirSampling_impl<T>::ReservoirSampling_impl(
	unsigned long sampleSize,
	Tools::Random& r, 
	unsigned long threshold
) : m_sampleSize(sampleSize / sizeof(unsigned long)),
    m_t(0),
    m_T(threshold),
    m_pRandom(&r)
{
}

template<class T> Sketches::ReservoirSampling_impl<T>::ReservoirSampling_impl(
	const byte* data,
	Tools::Random& r)
{
}

template<class T> Sketches::ReservoirSampling_impl<T>::~ReservoirSampling_impl()
{
}

template<class T> void Sketches::ReservoirSampling_impl<T>::insert(T key)
{
	static double W =
		std::exp(
			-std::log(m_pRandom->nextUniformDouble()) /
			static_cast<double>(m_sampleSize)
		);

	static bool bInsert = false;
	static unsigned long G = 0;

	if (bInsert)
	{
		replace(key);
		bInsert = false;
	}
	else if (m_t < m_sampleSize)
	{
		m_reservoir.push_back(key);
	}
	else if (G > 0)
	{
		G--;
		if (G == 0) bInsert = true;
	}
	else if (m_t <= m_T * m_sampleSize)
	{
		// Algorithm X
		double V = m_pRandom->nextUniformDouble();
		double t = m_t + 1;
		double quot =
			static_cast<double>(t) / static_cast<double>(t - m_sampleSize);

		while (quot > V)
		{
			G++;
			t++;
			quot *=
				static_cast<double>(t - m_sampleSize) / static_cast<double>(t);
		}

		if (G == 0) replace(key);
		else G--;
	}
	else
	{
		// Algorithm Z
		while (true)
		{
			unsigned long term = m_t - m_sampleSize + 1;
			double U = m_pRandom->nextUniformDouble();
			double X = static_cast<double>(m_t) * (W - 1.0);
			G = static_cast<unsigned long>(std::floor(X));
			double lhs =
				std::exp(std::log(((U * std::pow((m_t + 1.0) /
				static_cast<double>(term), 2.0)) * (term + G)) /
				(m_t + X)) / static_cast<double>(m_sampleSize));
			double rhs =
				(((m_t + X) / (static_cast<double>(term + G))) * term) /
				static_cast<double>(m_t);

			if (lhs < rhs)
			{
				W = rhs / lhs;
				break;
			}
		
			double y =
				(((U * (m_t + 1)) /	static_cast<double>(term)) *
				(m_t + G + 1.0)) / (m_t + X);
			unsigned long denom, numer_lim;

			if (m_sampleSize < G)
			{
				denom = m_t;
				numer_lim = term + G;
			}
			else
			{
				denom = m_t - m_sampleSize + G;
				numer_lim = m_t + 1;
			}
			
			for (unsigned long numer = m_t + G; numer >= numer_lim; numer--)
			{
				y = (y * numer) / static_cast<double>(denom);
				denom--;
			}
	
			W =
				std::exp(-std::log(m_pRandom->nextUniformDouble()) /
				static_cast<double>(m_sampleSize));

			if (
				std::exp(std::log(y) / static_cast<double>(m_sampleSize)) <=
				(m_t + X) / static_cast<double>(m_t)
			) break;
		}

		if (G == 0) replace(key);
		else G--;
	}

	m_t++;
}

template<class T> void Sketches::ReservoirSampling_impl<T>::erase(T key)
{
	throw Tools::NotSupportedException(
		"ReservoirSampling: Reservoir sampling does not support deletions."
	);
}

template<class T> void Sketches::ReservoirSampling_impl<T>::clear()
{
	m_t = 0;
	m_reservoir.clear();
}

template<class T> unsigned long
Sketches::ReservoirSampling_impl<T>::getInputLength() const
{
	return m_t;
}

template<class T> unsigned long
Sketches::ReservoirSampling_impl<T>::getFrequency(
	T id
) const
{
	double f = static_cast<double>(m_t) / static_cast<double>(m_sampleSize);
	
	unsigned long c = 0;
	for (unsigned long i = 0; i < m_reservoir.size(); i++)
	{
		if (m_reservoir[i] == id) c++;
	}

	return static_cast<unsigned long>(c * f);
}

template<class T> void Sketches::ReservoirSampling_impl<T>::replace(T key)
{
	unsigned long M = m_pRandom->nextUniformLong(0L, m_sampleSize);
	m_reservoir[M] = key;
}

/**************************************************/
/*            Template specializations            */
/**************************************************/

// This is the actual template. It does not have to
// be visible, since it cannot be instantiated in
// any case.
template<class T> class Sketches::ReservoirSampling
 : public Sketches::ReservoirSampling_impl<T>
{
	ReservoirSampling(
		double epsilon,
		double delta,
		Tools::Random& r,
		unsigned long threshold = 22
	);

	ReservoirSampling(
		unsigned long reservoirSize,
		Tools::Random& r,
		unsigned long threshold = 22
	);

	ReservoirSampling(const byte* data, Tools::Random& r);
	virtual ~ReservoirSampling();

	virtual unsigned long getSize() const = 0;
	virtual void getData(byte** data, unsigned long& length) const = 0;
};

Sketches::ReservoirSampling<unsigned long>::ReservoirSampling(
	double epsilon,
	double delta,
	Tools::Random& r,
	unsigned long threshold
) : ReservoirSampling_impl<unsigned long>(epsilon, delta, r, threshold)
{
}

Sketches::ReservoirSampling<unsigned long>::ReservoirSampling(
	unsigned long sampleSize,
	Tools::Random& r, 
	unsigned long threshold
) : ReservoirSampling_impl<unsigned long>(sampleSize, r, threshold)
{
}

Sketches::ReservoirSampling<unsigned long>::ReservoirSampling(
	const byte* data,
	Tools::Random& r
) : ReservoirSampling_impl<unsigned long>(data, r)
{
	memcpy(&m_sampleSize, data, sizeof(unsigned long));
	data += sizeof(unsigned long);
	memcpy(&m_t, data, sizeof(unsigned long));
	data += sizeof(unsigned long);
	memcpy(&m_T, data, sizeof(unsigned long));
	data += sizeof(unsigned long);

	unsigned long size;
	memcpy(&size, data, sizeof(unsigned long));
	data += sizeof(unsigned long);

	unsigned long l;
	for (unsigned long i = 0; i < size; i++)
	{
		memcpy(&l, data, sizeof(unsigned long));
		data += sizeof(unsigned long);
		m_reservoir.push_back(l);
	}

	m_pRandom = &r;
}

Sketches::ReservoirSampling<unsigned long>::~ReservoirSampling()
{
}

unsigned long Sketches::ReservoirSampling<unsigned long>::getSize() const
{
	return
		4 * sizeof(unsigned long) +
		m_reservoir.size() * sizeof(unsigned long);
}

void Sketches::ReservoirSampling<unsigned long>::getData(
	byte** data,
	unsigned long& length
) const
{
	length = getSize();
	*data = new byte[length];
	byte* p = *data;

	memcpy(p, &m_sampleSize, sizeof(unsigned long));
	p += sizeof(unsigned long);
	memcpy(p, &m_t, sizeof(unsigned long));
	p += sizeof(unsigned long);
	memcpy(p, &m_T, sizeof(unsigned long));
	p += sizeof(unsigned long);

	unsigned long l = m_reservoir.size();
	memcpy(p, &l, sizeof(unsigned long));
	p += sizeof(unsigned long);

	for (unsigned long i = 0; i < m_reservoir.size(); i++)
	{
		l = m_reservoir[i];
		memcpy(p, &l, sizeof(unsigned long));
		p += sizeof(unsigned long);
	}

	assert(p == (*data) + length);
}

Sketches::ReservoirSampling<std::string>::ReservoirSampling(
	double epsilon,
	double delta,
	Tools::Random& r,
	unsigned long threshold
) : ReservoirSampling_impl<std::string>(epsilon, delta, r, threshold)
{
}

Sketches::ReservoirSampling<std::string>::ReservoirSampling(
	unsigned long sampleSize,
	Tools::Random& r, 
	unsigned long threshold
) : ReservoirSampling_impl<std::string>(sampleSize, r, threshold)
{
}

Sketches::ReservoirSampling<std::string>::ReservoirSampling(
	const byte* data,
	Tools::Random& r
) : ReservoirSampling_impl<std::string>(data, r)
{
	memcpy(&m_sampleSize, data, sizeof(unsigned long));
	data += sizeof(unsigned long);
	memcpy(&m_t, data, sizeof(unsigned long));
	data += sizeof(unsigned long);
	memcpy(&m_T, data, sizeof(unsigned long));
	data += sizeof(unsigned long);

	unsigned long size;
	memcpy(&size, data, sizeof(unsigned long));
	data += sizeof(unsigned long);

	for (unsigned long i = 0; i < size; i++)
	{
		std::string s(reinterpret_cast<const char*>(data));
		data += s.size() + 1;
		m_reservoir.push_back(s);
	}

	m_pRandom = &r;
}

Sketches::ReservoirSampling<std::string>::~ReservoirSampling()
{
}

unsigned long Sketches::ReservoirSampling<std::string>::getSize() const
{
	unsigned long ret = 4 * sizeof(unsigned long);

	for (unsigned long i = 0; i < m_reservoir.size(); i++)
		ret += m_reservoir[i].size() + 1;

	return ret;
}

void Sketches::ReservoirSampling<std::string>::getData(
	byte** data,
	unsigned long& length
) const
{
	length = getSize();
	*data = new byte[length];
	byte* p = *data;

	memcpy(p, &m_sampleSize, sizeof(unsigned long));
	p += sizeof(unsigned long);
	memcpy(p, &m_t, sizeof(unsigned long));
	p += sizeof(unsigned long);
	memcpy(p, &m_T, sizeof(unsigned long));
	p += sizeof(unsigned long);

	unsigned long l = m_reservoir.size();
	memcpy(p, &l, sizeof(unsigned long));
	p += sizeof(unsigned long);

	for (unsigned long i = 0; i < m_reservoir.size(); i++)
	{
		memcpy(p, m_reservoir[i].c_str(), m_reservoir[i].size());
		p += m_reservoir[i].size();
		*p = 0;
		p++;
	}

	assert(p == (*data) + length);
}

