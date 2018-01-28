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

#ifndef __sketches_h
#define __sketches_h

#include <tools/Tools.h>

namespace Sketches
{
	enum HashType
	{
		HT_UNIVERSAL = 0x0,
		HT_SHA1
	};

	/**************************************************/
	/*                  LossyCounting                 */
	/**************************************************/
	// This implements the Lossy Counting algorithm presented in:
	// G. S. Manku, R. Motwani
	// Approximate Frequency Counts over Data Streams, VLDB 2002
	class LossyCounting
	{
	public:
		LossyCounting(double epsilon);
		virtual ~LossyCounting();

		virtual void insert(const std::string& id, unsigned long val = 1);
		virtual void erase(const std::string& id, unsigned long val = 1);
		virtual void insert(unsigned long id, unsigned long val = 1);
		virtual void erase(unsigned long id, unsigned long val = 1);
		virtual void clear();

		virtual unsigned long getNumberOfEntries() const;
			// returns the total number of entries in the frequent list.

		virtual unsigned long getMaxFrequency(const std::string& id) const;
			// returns an upper bound of the true frequency of element id.

		virtual unsigned long getMinFrequency(const std::string& id) const;
			// returns a lower bound of the true frequency of element id.

		virtual std::map<std::string, std::pair<unsigned long, unsigned long> >
		getFrequent(
			unsigned long theta
		) const;
			// returns all elements with true frequencies
			// PROBABLY higher than theta.
			// The key is the element identifier (id).
			// The data are the estimated frequency of the
			// element (ef) and the maximum possible error (delta),
			// such that for the true frequency
			// (f) it holds that: ef <= f <= ef + delta.

		virtual unsigned long getInputLength() const;
			// returns the total number of insertions
			// minus the total number of deletions from the sketch.

		virtual double getEpsilon() const;
			// returns the user specified accuracy parameter.

		virtual unsigned long getSize() const;
			// returns the total size of the sketch in bytes.

	private:
		class Entry
		{
		public:
			Entry(unsigned long c, unsigned long d);

			unsigned long m_count;
			unsigned long m_delta;
		};

		double m_epsilon;
		unsigned long m_w;
		unsigned long m_L;
		unsigned long m_bucket;
		std::map<std::string, Entry> m_sketch;
	};

	/**************************************************/
	/*                   Bloom Filter                 */
	/**************************************************/
	// A Bloom filter uses only one vector of
	// counters but hashes each item multiple times.
	// When SHA1 is used the domain is infinite
	// (SHA1 can hash strings of arbitrary length).
	// Nevertheless, only up to 2^16-1 counters and
	// at most 10 hash functions can be used.
	// When a Universal hash is used the domain is
	// sizeof(UniversalHash::value_type), but an arbitrary number
	// of counters and hash functions can be used.
	// The ids are converted to value_type using atoll.
	class BloomFilter
	{
	public:
		BloomFilter(
			unsigned long bits,
			unsigned long hashes,
			HashType t = Sketches::HT_UNIVERSAL
		);
		BloomFilter(
			unsigned long bits,
			unsigned long hashes,
			Tools::Random& r
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The random number generator is used to
			// produce the coefficients of the hash functions.

		BloomFilter(
			unsigned long bits,
			const std::vector<Tools::UniversalHash>& hashes
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The hash functions are explicitly provided.

		BloomFilter(const byte* data);
		virtual ~BloomFilter();

		virtual void insert(const std::string& id);
		virtual void erase(const std::string& id);
		virtual void insert(const Tools::UniversalHash::value_type& id);
		virtual void erase(const Tools::UniversalHash::value_type& id);
		virtual void clear();
	
		virtual bool contains(const std::string& id) const;
		virtual bool contains(const Tools::UniversalHash::value_type& id) const;

		virtual void merge(const BloomFilter& in);
		virtual BloomFilter getMerged(const BloomFilter& in) const;

		virtual void intersect(const BloomFilter& in);
		virtual BloomFilter getIntersection(const BloomFilter& in) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;
		virtual unsigned long getNumberOfBitsSet() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		HashType m_type;
		unsigned long m_bits;
		unsigned long m_hashes;
		std::vector<bool> m_filter;
		std::vector<Tools::UniversalHash> m_hash;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os,
			const Sketches::BloomFilter& s
		);
	};

	std::ostream& operator<<(
		std::ostream& os, const Sketches::BloomFilter& s
	);

	/**************************************************/
	/*            Multi-stage Bloom Filter            */
	/**************************************************/
	// A Multi-stage Bloom filter uses multiple
	// vectors of counters and a different
	// hash function per vector.
	// When SHA1 is used the domain is infinite
	// (SHA1 can hash strings of arbitrary length).
	// Nevertheless, only up to 2^16-1 counters and
	// at most 10 hash functions can be used.
	// When a Universal hash is used the domain is
	// sizeof(UniversalHash::value_type), but an arbitrary number
	// of counters and hash functions can be used.
	// The ids are converted to value_type using atoll.
	class MultiBloomFilter
	{
	public:
		MultiBloomFilter(
			unsigned long bits,
			unsigned long hashes,
			HashType t = HT_UNIVERSAL
		);
		MultiBloomFilter(
			unsigned long bits,
			unsigned long hashes,
			Tools::Random& r
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The random number generator is used to
			// produce the coefficients of the hash functions.

		MultiBloomFilter(
			unsigned long bits,
			const std::vector<Tools::UniversalHash>& hashes
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The hash functions are explicitly provided.

		MultiBloomFilter(const byte* data);
		virtual ~MultiBloomFilter();

		virtual void insert(const std::string& id);
		virtual void erase(const std::string& id);
		virtual void insert(const Tools::UniversalHash::value_type& id);
		virtual void erase(const Tools::UniversalHash::value_type& id);
		virtual void clear();
	
		virtual bool contains(const std::string& id) const;
		virtual bool contains(const Tools::UniversalHash::value_type& id) const;

		virtual void merge(const MultiBloomFilter& in);
		virtual MultiBloomFilter getMerged(
			const MultiBloomFilter& in
		) const;

		virtual void intersect(const MultiBloomFilter& in);
		virtual MultiBloomFilter getIntersection(
			const MultiBloomFilter& in
		) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;
		virtual unsigned long getNumberOfBitsSet() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		HashType m_type;
		unsigned long m_bits;
		unsigned long m_hashes;
		std::vector<bool> m_filter;
		std::vector<Tools::UniversalHash> m_hash;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os,
			const Sketches::MultiBloomFilter& s
		);
	};

	std::ostream& operator<<(
		std::ostream& os, const Sketches::MultiBloomFilter& s
	);

	/**************************************************/
	/*              Counting Bloom Filter             */
	/**************************************************/
	// This is an implementation of a counting Bloom filter from:
	// L. Fan, P. Cao, J. Almeida and A. Z. Broder.
	// Summary Cache: A Scalable Wide-Area Web Cache Sharing Protocol.
	// IEEE/ACM Transactions on Networking 8(3), 2000.
	//
	// Every vector entry is a 4-bit counter instead
	// of a bit. Counters do not overflow or underflow.
	// A check is performed.
	// When SHA1 is used the domain is infinite
	// (SHA1 can hash strings of arbitrary length).
	// Nevertheless, only up to 2^16-1 counters and
	// at most 10 hash functions can be used.
	// When a Universal hash is used the domain is
	// sizeof(UniversalHash::value_type), but an arbitrary number
	// of counters and hash functions can be used.
	// The ids are converted to value_type using atoll.
	//
	// NOTICE: The usefulness of counting Bloom filters
	// with counters larger than 4 bits is "iffy".
	// For dense sets, it is possible that in order to
	// create an accurate filter with larger counters,
	// you need as much space as an exact solution would use
	// (see A. Broder and M. Mitzenmacher.
	// Network Applications of Bloom Filters: A Survey).
	// NOTICE: If you would like to use a similar sketch
	// with larger counters, please use the FastAMS sketch.
	class CountingBloomFilter
	{
	public:
		CountingBloomFilter(
			unsigned long counters,
			unsigned long hashes,
			HashType t = HT_UNIVERSAL
		);
		CountingBloomFilter(
			unsigned long counters,
			unsigned long hashes,
			Tools::Random& r
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The random number generator is used to
			// produce the coefficients of the hash functions.

		CountingBloomFilter(
			unsigned long counters,
			const std::vector<Tools::UniversalHash>& hashes
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The hash functions are explicitly provided.

		CountingBloomFilter(const CountingBloomFilter& in);
		CountingBloomFilter(const byte* data);
		virtual ~CountingBloomFilter();

		virtual CountingBloomFilter& operator=(const CountingBloomFilter& in);

		virtual void insert(const std::string& id, byte val = 1);
		virtual void erase(const std::string& id, byte val = 1);
		virtual void insert(
			const Tools::UniversalHash::value_type& id,
			byte val = 1
		);
		virtual void erase(
			const Tools::UniversalHash::value_type& id,
			byte val = 1
		);
		virtual void clear();
	
		virtual byte getFrequency(const std::string& id) const;
		virtual byte getFrequency(
			const Tools::UniversalHash::value_type& id
		) const;
	
		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		HashType m_type;
		unsigned long m_counters;
		unsigned long m_hashes;
		byte* m_pFilter;
		std::vector<Tools::UniversalHash> m_hash;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os,
			const Sketches::CountingBloomFilter& s
		);
	};

	std::ostream& operator<<(
		std::ostream& os,
		const Sketches::CountingBloomFilter& s
	);

	/**************************************************/
	/*       Multi-stage Counting Bloom Filter        */
	/**************************************************/
	// This is an implementation of a multi-stage counting
	// Bloom filter, where every vector entry is a 4-bit
	// counter instead of a bit. Counters do not overflow
	// or underflow. A check is performed.
	// When SHA1 is used the domain is infinite
	// (SHA1 can hash strings of arbitrary length).
	// Nevertheless, only up to 2^16-1 counters and
	// at most 10 hash functions can be used.
	// When a Universal hash is used the domain is
	// sizeof(UniversalHash::value_type, but an arbitrary number
	// of counters and hash functions can be used.
	// The ids are converted to value_type using atoll.
	//
	// NOTICE: The usefulness of counting Bloom filters
	// with counters larger than 4 bits is "iffy".
	// For dense sets, it is possible that in order to
	// create an accurate filter with larger counters,
	// you need as much space as an exact solution would use
	// (see A. Broder and M. Mitzenmacher.
	// Network Applications of Bloom Filters: A Survey).
	// NOTICE: If you would like to use a similar sketch
	// with larger counters, please use the FastAMS sketch.
	class MultiCountingBloomFilter
	{
	public:
		MultiCountingBloomFilter(
			unsigned long counters,
			unsigned long hashes,
			HashType t = HT_UNIVERSAL
		);
		MultiCountingBloomFilter(
			unsigned long counters,
			unsigned long hashes,
			Tools::Random& r
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The random number generator is used to
			// produce the coefficients of the hash functions.

		MultiCountingBloomFilter(
			unsigned long counters,
			const std::vector<Tools::UniversalHash>& hashes
		);
			// When using this constructor the hash type
			// of the filter is always HT_UNIVERSAL.
			// The hash functions are explicitly provided.

		MultiCountingBloomFilter(const MultiCountingBloomFilter& in);
		MultiCountingBloomFilter(const byte* data);
		virtual ~MultiCountingBloomFilter();

		virtual MultiCountingBloomFilter& operator=(
			const MultiCountingBloomFilter& in
		);

		virtual void insert(const std::string& id, byte val = 1);
		virtual void erase(const std::string& id, byte val = 1);
		virtual void insert(
			const Tools::UniversalHash::value_type& id,
			byte val = 1
		);
		virtual void erase(
			const Tools::UniversalHash::value_type& id,
			byte val = 1
		);
		virtual void clear();
	
		virtual byte getFrequency(const std::string& id) const;
		virtual byte getFrequency(
			const Tools::UniversalHash::value_type& id
		) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		HashType m_type;
		unsigned long m_counters;
		unsigned long m_hashes;
		byte* m_pFilter;
		unsigned long m_filterSize;
		std::vector<Tools::UniversalHash> m_hash;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os,
			const MultiCountingBloomFilter& s
		);
		friend class CountMin;
	};

	std::ostream& operator<<(
		std::ostream& os,
		const Sketches::MultiCountingBloomFilter& s
	);

	/**************************************************/
	/*                 CountMin Sketch                */
	/**************************************************/
	// This implements the CountMin sketch as proposed in:
	// R. Motwani and P. Raghavan.
	// Randomized Algorithms
	// Cambridge International Series on Parallel Computation.
	// and
	// G. Cormode and S. Muthukrishnan
	// An Improved Data Stream Summary: The Count-Min Sketch and
	// its Applications, Journal of Algorithms 55(1), 2005
	//
	// The CountMin sketch is the same as a Multi-stage
	// Bloom filter. It provides guarantees on the returned
	// answers by enforcing specific vector sizes and number
	// of hash functions.
	// NOTICE: If you would like to use a similar sketch
	// with larger counters, please use the FastAMS sketch.
	class CountMin : public MultiCountingBloomFilter
	{
	public:
		virtual ~CountMin();

		CountMin(double epsilon, double delta);
		CountMin(double epsilon, double delta, Tools::Random& r);
	};

	/**************************************************/
	/*           FM and Summation FM Sketch           */
	/**************************************************/
	// This implements the PCSA FM sketch as described in:
	// P. Flajolet and G. N. Martin.
	// Probabilistic Counting Algorithms for data base applications.
	// Journal of Computer and System Sciences, 31, 1985
	//
	// It also implements the sumation FM sketch which
	// supports efficient execution of "bulk insertions" to
	// the FM sketch, where each insertion has weight larger
	// than one. The sumation FM code has been kindly contributed
	// by the authors of:
	// J. Considine, F. Li, G. Kollios, J. Byers
	// Approximate aggregation techniques for sensor databases, ICDE 2004
	//
	// For HT_UNIVERSAL the item ids are converted
	// from std::string to UniversalHash::value_type using atoll.
	// For HT_SHA1 the domain size is infinite,
	// but the HT_SHA1 version does not support bulk insertions.
	// It inserts the values one by one iteratively.
	// Below is a table of supported domain sizes and sketch
	// sizes, depending on the hash function used:
	//
	// Hash type      Domain size     Bits per Bitmap    Bitmaps
	// HT_UNIVERSAL   value_type      infinite           infinite
	// HT_SHA1        infinite        152                256
	class FM
	{
	public:
		FM(
			unsigned long bits = 32,
			unsigned long bitmaps = 64,
			HashType type = HT_SHA1
		);
		FM(const byte* data);
		virtual ~FM();

		virtual void insert(const std::string& id, unsigned long val = 1);
		virtual void erase(const std::string& id, unsigned long val = 1);
		virtual void insert(
			const Tools::UniversalHash::value_type& id,
			unsigned long val = 1
		);
		virtual void erase(
			const Tools::UniversalHash::value_type& id,
			unsigned long val = 1
		);
		virtual void clear();

		virtual unsigned long getCount() const;

		virtual void merge(const FM& f);
			// merge this instance with the input FM and store in this instance.
		virtual void merge(const FM& f1, const FM& f2);
			// merge the input FMs and store the result in this instance.
		virtual FM getMerged(const FM& f) const;
			// merge this instance with the input FM
			// and store the result only on the returned value.
		virtual void reset(const byte* data);
		virtual bool isSubsumedBy(const FM& f) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfVectors() const;
		virtual unsigned long getUncompressedSize() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		unsigned long pickOffset(Tools::Random& r) const;
		unsigned long fillLength(unsigned long val) const;
		long pickBinomial(long nIn, double pIn, Tools::Random& rand) const;

		static const double PHI = 0.77351;

		HashType m_type;
		std::vector<std::vector<bool> > m_bitmap;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os, const Sketches::FM& s
		);
	};

	std::ostream& operator<<(std::ostream& os, const Sketches::FM& s);

	/**************************************************/
	/*                    AMS Sketch                  */
	/**************************************************/
	// The AMS sketch as it appeared in:
	// N. Alon, Y. Matias, and M. Szegedy.
	// The space complexity of approximating the frequency moments. STOC 96.
	//
	// Using SHA1 here does not make sense since we
	// need fourwise independent hashing.
	// Ideally, we should be storing only the PRG seed
	// and produce all neccessary hashes on the fly here,
	// but I will store the hashes for efficiency and
	// conviniency.
	class AMS
	{
	public:
		AMS(double e, double d);
		AMS(unsigned long s1, unsigned long s2);
		AMS(double e, double d, unsigned long seed);
		AMS(unsigned long s1, unsigned long s2, unsigned long seed);
		AMS(const AMS& in);
		AMS(const byte* data);
		virtual ~AMS();

		virtual AMS& operator=(const AMS& in);

		virtual void insert(const std::string& id, long val = 1);
		virtual void erase(const std::string& id, long val = 1);
		virtual void insert(
			const Tools::UniversalHash::value_type& id,
			long val = 1
		);
		virtual void erase(
			const Tools::UniversalHash::value_type& id,
			long val = 1
		);
		virtual void clear();

		virtual double getF2Norm() const;
		virtual long getFrequency(const std::string& id) const;
		virtual long getFrequency(
			const Tools::UniversalHash::value_type& id
		) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfVectors() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& len) const;

	public:
		unsigned long m_seed;
		unsigned long m_s1;
		unsigned long m_s2;
		long* m_pAtomicSketch;
		std::vector<Tools::UniversalHash> m_hash;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os, const Sketches::AMS& s
		);
	};

	std::ostream& operator<<(std::ostream& os, const Sketches::AMS& s);

	/**************************************************/
	/*                Fast AMS Sketch                 */
	/**************************************************/
	// The FastAMS sketch as proposed in:
	// Graham Cormode and Minos Garofalakis.
	// Sketching Streams Through the Net: Distributed
	// Approximate Query Tracking. VLDB 2005.
	//
	// Using SHA1 here does not make sense since we
	// need fourwise independent hashing.
	// Ideally, we should be storing only the PRG seed
	// and produce all neccessary hashes on the fly here,
	// but I will store the hashes for efficiency and
	// conviniency.
	class FastAMS
	{
	public:
		FastAMS(unsigned long counters, unsigned long hashes);
		FastAMS(
			unsigned long counters,
			unsigned long hashes,
			unsigned long seed
		);
		FastAMS(const FastAMS& in);
		FastAMS(const byte* data);
		virtual ~FastAMS();

		virtual FastAMS& operator=(const FastAMS& in);

		virtual void insert(const std::string& id, long val = 1);
		virtual void erase(const std::string& id, long val = 1);
		virtual void insert(
			const Tools::UniversalHash::value_type& id,
			long val = 1
		);
		virtual void erase(
			const Tools::UniversalHash::value_type& id,
			long val = 1
		);
		virtual void clear();

		virtual long getFrequency(const std::string& id) const;
		virtual long getFrequency(
			const Tools::UniversalHash::value_type& id
		) const;

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& len) const;

	private:
		unsigned long m_seed;
		unsigned long m_counters;
		long* m_pFilter;
		std::vector<Tools::UniversalHash> m_hash;
			// this is used to hash a specific element
			// on a specific counter per row.
		std::vector<Tools::UniversalHash> m_fourwiseHash;
			// this is used to hash a specific element
			// to -1 or +1 per row.

		friend std::ostream& Sketches::operator<<(
				std::ostream& os, const Sketches::FastAMS& s
		);
	};

	std::ostream& operator<<(std::ostream& os, const Sketches::FastAMS& s);

	/**************************************************/
	/*                Reservoir Sampling              */
	/**************************************************/
	// This implements the reservoir sampling algorithm given in
	// J. S. Vitter
	// Random Sampling with a Reservoir, TOMS, Vol. 11, Issue 1, 1985
	//
	// I will put all the funcitonality in ReservoirSampling_impl
	// so that I don't have to duplicate all the code for
	// every template specialization.
	// Actually, there is no need to use template specializations
	// here. Direct subclasses of ReservoirSampling_impl
	// would do just fine. I use specializations just
	// because they look better when declaring them.
	template<class T> class ReservoirSampling;

	template<class T> class ReservoirSampling_impl
	{
	public:
		ReservoirSampling_impl(
			double epsilon,
			double delta,
			Tools::Random& r,
			unsigned long threshold = 22
		);
			// the sample size will be determined according
			// to probabilistic guarantees.

		ReservoirSampling_impl(
			unsigned long reservoirSize,
			Tools::Random& r,
			unsigned long threshold = 22
		);
			// the reservoir size is given by the caller in number of bytes.

		ReservoirSampling_impl(const byte* data, Tools::Random& r);
		virtual ~ReservoirSampling_impl();

		virtual void insert(T id);
		virtual void erase(T id);
		virtual void clear();

		virtual unsigned long getInputLength() const;
		virtual unsigned long getFrequency(T id) const;

		virtual unsigned long getSize() const = 0;
		virtual void getData(byte** data, unsigned long& length) const = 0;

	private:
		void replace(T key);

		unsigned long m_sampleSize;
		unsigned long m_t;
		unsigned long m_T;
		std::vector<T> m_reservoir;
		Tools::Random* m_pRandom;
			// Assume that a system wide source of
			// randomness is available in general.

		friend class ReservoirSampling<unsigned long>;
		friend class ReservoirSampling<std::string>;
	};

	template<> class ReservoirSampling<unsigned long>
	 : public ReservoirSampling_impl<unsigned long>
	{
	public:
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

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;
	};

	template<> class ReservoirSampling<std::string>
	 : public ReservoirSampling_impl<std::string>
	{
	public:
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

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;
	};

	/**************************************************/
	/*                 Counting Sample                */
	/**************************************************/
	// This implements the counting sampling algorithm of:
	// P. B. Gibbons and Y. Matias
	// New sampling-based summary statistics for improving
	// approximate query answers, SIGMOD, 1998
	template<class T> class CountingSample;

	template<class T> class CountingSample_impl
	{
	public:
		CountingSample_impl(unsigned long reservoirSize, Tools::Random& r);
			// the reservoir size is given by the caller
			// in number of entries in the sample.
		CountingSample_impl(const CountingSample_impl& in);
		CountingSample_impl(const byte* data, Tools::Random& r);
		virtual ~CountingSample_impl();

		virtual CountingSample_impl& operator=(const CountingSample_impl& in);

		virtual void insert(T key);
		virtual void erase(T key);
		virtual void clear();

		virtual unsigned long getInputLength() const;
		virtual double getSamplingFactor() const;
		virtual unsigned long getFrequency(T id) const;
		virtual std::map<T, unsigned long> getFrequent(
			unsigned long theta
		) const;
		virtual std::vector<std::pair<T, unsigned long> > getSample() const;
		virtual unsigned long getTotalKeyTallies() const;

		virtual unsigned long getSize() const = 0;
		virtual void getData(byte** data, unsigned long& length) const = 0;

	private:
		void resample();

		class Key
		{
		public:
			Key(T key) : m_key(key) {}
			virtual ~Key() {}
			T m_key;
		};

		class KeyTally : public Key
		{
		public:
			KeyTally(T key, unsigned long c) : Key(key), m_count(c) {}
			unsigned long m_count;
		};

		struct KeyLess : public std::binary_function<const Key*, const Key*, bool>
		{
			bool operator()(const Key* k1, const Key* k2) const
			{
				return k1->m_key < k2->m_key;
			}
		};

		unsigned long m_N;
		unsigned long m_sampleSize;
		double m_tau;
		std::set<Key*, KeyLess> m_reservoir;
		Tools::Random* m_pRandom;
			// the Mersenne generator needs to store a 19968
			// bit string as a seed, which might be an overkill
			// compared to the sample size in some cases.
			// It is much better to assume that a system wide
			// source of randomness is available in general.

		friend class CountingSample<unsigned long>;
		friend class CountingSample<std::string>;
	};

	template<> class CountingSample<unsigned long>
	 : public CountingSample_impl<unsigned long>
	{
	public:
		CountingSample(unsigned long reservoirSize, Tools::Random& r);
		CountingSample(const byte* data, Tools::Random& r);
		CountingSample(const CountingSample<unsigned long>& in);
		virtual ~CountingSample();
		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;
	};

	template<> class CountingSample<std::string>
	 : public CountingSample_impl<std::string>
	{
	public:
		CountingSample(unsigned long reservoirSize, Tools::Random& r);
		CountingSample(const byte* data, Tools::Random& r);
		CountingSample(const CountingSample<std::string>& in);
		virtual ~CountingSample();
		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;
	};

	/**************************************************/
	/*                 Quantile Digest                */
	/**************************************************/
	// This implements the quantile digest algorithm given in:
	// S. Nath and P. B. Gibbons and S. Seshan and Z. R. Anderson
	// Synopsis diffusion for robust aggregation in sensor networks,
	// SENSYS 2004
	//
	// TODO: implement insert.
	class QuantileDigest
	{
	public:
		QuantileDigest(
			unsigned long k,
			long min,
			long max,
			const std::vector<long>& data
		);
			// k is the compression parameter.
			// min is the lower bound of the value domain
			// max is the upper bound of the value domain
			// (these should be the same across all Quantile
			// Digests to be merged)
			// data is the input array of values.

		QuantileDigest(
			unsigned long k,
			long min,
			long max,
			const std::map<long, unsigned long>& data
		);
				// data is a histogram of the input array of values

		QuantileDigest(const QuantileDigest& in);
		QuantileDigest(const byte* data);
		virtual ~QuantileDigest();

		virtual QuantileDigest& operator=(const QuantileDigest& in);

		virtual void getQuantile(
			double q,
			long& value,
			unsigned long& rank
		) const;

		virtual unsigned long getInputLength() const;

		virtual void merge(const QuantileDigest& in);
		virtual QuantileDigest getMerged(const QuantileDigest& in) const;

		virtual unsigned long getSize() const;
		virtual void getData(byte** data, unsigned long& length) const;

	private:
		class Node
		{
		public:
			Node(long low, long high, unsigned long c);

			struct PostOrder : public std::binary_function<
				const Node*, const Node*, bool
			>
			{
				bool operator()(const Node* n1, const Node* n2);
			};

			// The rank of the node in the post order traversal can be used
			// here as a unique node id to save space. Nevertheless, the
			// computation is expensive and hence avoided.
			// It can be computed as follows:
			// while rank != 2^y - 1, for any y:
			//   x = floor(log_2(rank));
			//   rank = rank - 2^x
			//   c += 2^(x-1)
			// low = c + 1
			// ...
			unsigned long m_low, m_high, m_c;
		};

		void initialize(
			const std::map<long, unsigned long>& data
		) throw(Tools::IllegalArgumentException);

		void compress(std::vector<std::set<Node*, Node::PostOrder> >& qd);

		long m_min;
		long m_max;
		unsigned long m_k;
		unsigned long m_N;
		std::vector<std::set<Node*, Node::PostOrder> > m_qd;

		friend std::ostream& Sketches::operator<<(
			std::ostream& os,
			const Sketches::QuantileDigest& s
		);
		friend class QuantileDigestFM;
	};

	std::ostream& operator<<(
		std::ostream& os,
		const Sketches::QuantileDigest& s
	);

	/**************************************************/
	/*                CountMin FM Sketch              */
	/**************************************************/
	// TODO: implement getData
	class CountMinFM
	{
	public:
		CountMinFM(
			const std::string& id,
			const std::map<std::string, unsigned long>& m,
			unsigned long counters,
			unsigned long hashes,
			unsigned long fmSize,
			unsigned long fmBitmaps,
			Tools::Random& r
		);
		CountMinFM(
			const std::string& id,
			const std::map<std::string, unsigned long>& m,
			unsigned long counters,
			unsigned long hashes,
			unsigned long fmSize,
			unsigned long fmBitmaps,
			HashType t
		);
		CountMinFM(
			unsigned long id,
			const std::map<std::string, unsigned long>& m,
			unsigned long counters,
			unsigned long hashes,
			unsigned long fmSize,
			unsigned long fmBitmaps,
			Tools::Random& r
		);
		CountMinFM(
			unsigned long id,
			const std::map<std::string, unsigned long>& m,
			unsigned long counters,
			unsigned long hashes,
			unsigned long fmSize,
			unsigned long fmBitmaps,
			HashType t
		);
		virtual ~CountMinFM();

		virtual unsigned long getFrequency(const std::string& id) const;

		virtual void merge(const CountMinFM& in);

		virtual unsigned long getVectorLength() const;
		virtual unsigned long getNumberOfHashes() const;

		virtual unsigned long getSize() const;

	private:
		void initialize(
			const std::string& id,
			const std::map<std::string, unsigned long>& m,
			unsigned long counters,
			unsigned long hashes,
			unsigned long fmSize,
			unsigned long fmBitmaps,
			HashType t,
			Tools::Random& r
		);

		HashType m_type;
		std::string m_id;
			// Every CountMinFM should have a unique id,
			// for proper merging.
		unsigned long m_counters;
		unsigned long m_hashes;
		std::vector<FM> m_filter;
		std::vector<Tools::UniversalHash> m_hash;
	};

	/**************************************************/
	/*             Quantile Digest FM                 */
	/**************************************************/
	// TODO: implement getData
	class QuantileDigestFM
	{
	public:
		QuantileDigestFM(
			const std::string& id,
			const QuantileDigest& in,
			unsigned long fmSize = 15,
			unsigned long fmBitmaps = 20
		);
		QuantileDigestFM(
			unsigned long id,
			const QuantileDigest& in,
			unsigned long fmSize = 15,
			unsigned long fmBitmaps = 20
		);
		QuantileDigestFM(const QuantileDigestFM& in);
		virtual ~QuantileDigestFM();

		virtual QuantileDigestFM& operator=(const QuantileDigestFM& in);

		virtual void getQuantile(
			double q, long& value, unsigned long& rank
		) const;
		virtual void getQuantile(
			unsigned long total, double q, long& value, unsigned long& rank
		) const;

		virtual void merge(const QuantileDigestFM& in);
		virtual QuantileDigestFM getMerged(const QuantileDigestFM& in) const;

		virtual unsigned long getInputLength() const;

		virtual unsigned long getSize() const;

	private:
		class Node
		{
		public:
			Node(
				const std::string& id,
				long low,
				long high,
				unsigned long rank,
				unsigned long c,
				unsigned long fmSize = 15,
				unsigned long fmBitmaps = 20
			);
			Node(
				long low,
				long high,
				const Sketches::FM& c
			);
 
			struct PostOrder : public std::binary_function<
				const Node*,
				const Node*,
				bool
			>
			{
				bool operator()(const Node* n1, const Node* n2);
			};

			static unsigned long getNodeRank(
				unsigned long low, unsigned long high
			);

			// FIXME: normally, I should be able to figure
			// out low and high by using the unique identifier
			// assigned by the post-order traversal of the tree
			// in order to save some space.
			unsigned long m_low, m_high;
			Sketches::FM m_c;
		};

		void initialize(
			const QuantileDigest& in,
			unsigned long fmSize,
			unsigned long fmBitmaps
		);
		void compress(std::vector<std::set<Node*, Node::PostOrder> >& qd);

		std::string m_id;
			// Every QuantileDigestFM should have a unique id,
			// for proper merging.
		long m_min;
		long m_max;
		unsigned long m_k;
		Sketches::FM m_N;
		std::vector<std::set<Node*, Node::PostOrder> > m_qd;
	};

	/**************************************************/
	/*                 Useful functions               */
	/**************************************************/

	template<class T> T getMedian(std::multiset<T>& v)
	{
		typename std::multiset<T>::iterator it = v.begin();
		unsigned long r = static_cast<unsigned long>(std::ceil(static_cast<double>(v.size()) / 2.0));

		for (unsigned long i = 1; i < r; i++) it++;

		if (v.size() % 2 != 0)
		{
			return *it;
		}
		else
		{
			T a = *it; it++;
			T b = *it;
			return static_cast<T>((a + b) / 2.0);
		}
	}
};

#endif /* __sketches_h */

