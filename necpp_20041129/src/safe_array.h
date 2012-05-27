/*
	Copyright (C) 2004  Timothy C.A. Molteno
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __safe_array__
#define __safe_array__

#include <iostream>
#include <string>
#include <sstream>

#ifdef NEC_ERROR_CHECK

#include "nec_exception.h"

class BoundsViol : public nec_exception
{
public:
	BoundsViol(const char* message, long index, long bound)
		: nec_exception(message)
	{
		m_message << "array index: " << index << " exceeds " << bound << std::endl;
	}
};
#endif

/**
	A Safe Array class for nec2++ that performs bounds checking if
	the macro NEC_ERROR_CHECK is defined at compile time.
	
	This class also includes some utility functions for handling
	common vector operations.
*/
template<typename T>
class safe_array
{
public:
	safe_array()
		: len_(0), rows_(0), cols_(0), data_(NULL), data_size_(0), own_data_(true)
	{ }

	safe_array(long size)
		: len_(0), rows_(0), cols_(0), data_(NULL), data_size_(0), own_data_(true)
	{
		resize(size);
	}
	
	safe_array(long n_rows, long n_cols)
		: len_(0), rows_(0), cols_(0), data_(NULL), data_size_(0), own_data_(true)
	{
		resize(n_rows, n_cols);
	}
	
	safe_array(const safe_array<T>& in_array)
		: len_(0), rows_(0), cols_(0), data_(NULL), data_size_(0), own_data_(true)
	{
		copy(in_array);
	}
	
	~safe_array()
	{
		if (own_data_)
			delete[] data_;
	}
	
	long size() const
	{
		return len_;
	}

	void resize(long n_rows, long n_cols)
	{
		rows_ = n_rows;
		cols_ = n_cols;
		
		resize(rows_ * cols_);
	}

	// Copy the contents of in_array to our array
	// resizing as appropriate.
	void copy(const safe_array<T>& in_array)
	{
		if (in_array.rows_ == 0)
			resize(in_array.len_);
		else
			resize(in_array.rows_,in_array.cols_);
			
		for (long i=0; i<len_; i++)
			data_[i] = in_array[i];
	}

	void resize(long new_length)
	{
#ifdef NEC_ERROR_CHECK
		if (! own_data_)
			throw new nec_exception("attempt to resize data we do not own");
#endif			
		if (new_length > data_size_)
		{
			T* new_data_ = new T[new_length];
			data_size_ = new_length;	// keep the actual data size...
			
			for (long i=0; i<len_; i++)
				new_data_[i] = data_[i];
		
			delete[] data_;
			data_ = new_data_;
			len_ = new_length;
		}
		else
		{
			len_ = new_length;
		}
	}
	
	// return the largest element of the array
	T max()
	{
		T ret = data_[check(0)];
		
		for (long i = 1; i < len_; i++ )
		{
			if ( data_[check(i)] > ret)
				ret = data_[check(i)];
		}
		return ret;
	}

	// return the sum of all elements in the array
	T sum(long start_index, long stop_index)
	{
		T ret = data_[check(start_index)];
		
		for (long i = start_index+1; i < stop_index; i++ )
		{
			ret += data_[check(i)];
		}
		return ret;
	}

	// return the sum of all elements in the array
	T sum()
	{
		return sum(0,len_);
	}

	// fill all elements of the array with x
	void fill(long start, long N, const T& x)
	{
		long stop = start + N;
		for (long i = start; i < stop; i++ )
		{
			data_[check(i)] = x;
		}
	}
	
	// fill all elements of the array with x
	void fill(const T& x)
	{
		fill(0,len_,x);
	}
		
	T& get(long row, long col)
	{
		return data_[check(row,col)];
	}

	void set(long row, long col, const T& x)
	{
		data_[check(row,col)] = x;
	}
	
	const T& operator[](long i) const
	{
		return data_[check(i)];
	}
	
	T& operator[](long i)
	{
		return data_[check(i)];
	}
	
	// if end_index is -1, then finish at the end of the array
	safe_array<T> sub_array(long start_index, long end_index = -1)
	{
		if (-1 == end_index)
			end_index = len_;
			
		return safe_array<T>(*this, start_index, end_index, false);
	}
	
	
	T* get_ptr() const
	{
		return data_;
	}
	
	
	safe_array<T>& operator=(const safe_array<T>& in_array)
	{
		copy(in_array);
		return *this;
	}

private:
	long len_;
	long rows_;
	long cols_;
	
	T*  data_;
	long data_size_;
	
	bool own_data_;
	
	// Used to constructing sub_array's
	safe_array(safe_array<T>& in_array, long start_index, long end_index, bool in_copy_data)
	{
		len_ = (end_index - start_index);
		rows_ = 0;
		cols_ = 0;
		
		if (in_copy_data)
		{
			T* data_ = new T[len_];
			data_size_ = len_;
			
			for (long i=0; i<len_; i++)
				data_[check(i)] = in_array[start_index + i];
		
			own_data_ = true;
		}
		else
		{
			data_ = in_array.get_ptr() + start_index;
			data_size_ = 0;
			own_data_ = false;
		}
	}
	
	inline long check(long i) const
	{
#ifdef NEC_ERROR_CHECK
		if (i < 0 || i >= len_)
			throw new BoundsViol("safe_array", i, len_);
#endif
		return i;
	}

	inline long check(long row, long col) const
	{
#ifdef NEC_ERROR_CHECK
		if (row < 0 || row >= rows_)
			throw new BoundsViol("safe_array", row, rows_);
		if (col < 0 || col >= cols_)
			throw new BoundsViol("safe_array", col, cols_);
#endif
		return check(row*cols_ + col);
	}
}; 


#endif /* __safe_array__ */
