/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017-2018 Ingo Weyrich <heckflosse67@gmx.de>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
void findMinMaxPercentile(const float* data, size_t size, float minPrct, float* minOut, float maxPrct, float* maxOut, int threads);
/* Same as findMinMaxPercentile but skips the initial min/max scan when the
 * caller already knows the data's range. minVal/maxVal must bound every
 * sample; passing the actual data range gives the most accurate histogram. */
void findMinMaxPercentileFromMinMax(const float* data, size_t size, float minVal, float maxVal, float minPrct, float* minOut, float maxPrct, float* maxOut, int threads);
#ifdef __cplusplus
}
#endif
