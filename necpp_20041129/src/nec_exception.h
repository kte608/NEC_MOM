/***************************************************************************
 *   Copyright (C) 2004 by Tim Molteno,,,                                  *
 *   tim@cyberiad                                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef __nec_exception__
#define __nec_exception__

#include <string>
#include <ostream>

class nec_exception
{
public:
	nec_exception(const char* message)
	{
		m_message << message;
	}

	nec_exception(const char* message, int code)
	{
		m_message << message << code;
	}
	
	void append(const char* message)
	{
		m_message << message;
	}
	
	void append(const int i)
	{
		m_message << i;
	}
	
	std::string get_message()
	{
		std::string ret = m_message.str();
		return ret;
	}
	
protected:
	std::stringstream m_message;
};


#endif /* __nec_exception__ */
