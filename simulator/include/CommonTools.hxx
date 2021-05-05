/*
 * CommonTools.hxx
 *
 *  Created on: 18 May 2017
 *      Author: futyand
 */

#include <vector>
#include <cmath>
#include <algorithm>

#ifndef SIMULATOR_INCLUDE_COMMONTOOLS_HXX_
#define SIMULATOR_INCLUDE_COMMONTOOLS_HXX_

namespace CommonTools {

	/** *************************************************************************
	 *  @brief Returns the median value of a vector
	 *
	 *  @param [in] values  Vector of values for which the median is to be
	 *  					calculated
	 */
	template<typename T> double median(std::vector<T> & values) {

		unsigned size = values.size();
		if (size==0) return 0;
		std::sort(values.begin(),values.end());

		if (size%2 == 0) {
			return (values[(size/2)-1] + values[size/2])/2.;
		} else {
			return values[size/2];
		}

	}

}

#endif /* SIMULATOR_INCLUDE_COMMONTOOLS_HXX_ */
