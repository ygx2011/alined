# ALineD: Augmented Line-Based DLT

This library is based on the paper "Pose Estimation from Line Correspondence using Direct Linear Transform" by Pribyl,B. 2016,
the R_and_T algorithm as shown in "Robust Methods for Estimating Pose and a Sensitivity Analysis" by Kumar,R.,Hanson,A.R.,1994,
loss functions as explained in "Bundle Adjustment - A modern Synthesis", by Triggs,B., 1999,
and the paper "A Brief Description of the Levenberg-Marquardt Algorithm Implemened by levmar", by Manolis,I., Lourakis, A., 2005.
It consists of a fully line-based DLT algorithm as described in the above paper by Pribyl, a modification of the R_and_T algorithm
and an Extended Kalman Filter for sensor fusion.

![alt text](https://github.com/AndreaLampart/alined/blob/master/img/alined.png "Pose from lines")


## Motivation

Modern SLAM systems use point feature correspondences to extract camera motion. Points are limiting in a way that they
contain less structural information about the environment. Not many fully line-based systems have been published.
This library should be a help to whoever will takle the problem of line-based SLAM or VO.

## Installation

The software is self-contained and uses no additional dependencies except Eigen3. First, you should install the Eigen library by following their steps provided on the official website. Afterwards just clone the repository:

```bash
git clone directory-name.git
```


## How to use

The library can be included into any project and can be run from within a ros node. The following will calculate the camera pose:


```c++
// Configuration is sent via bitmasking
Alined alined(AL_COMBINED_LINES|AL_USE_REFINE|AL_LEVENBERG_MARQUARDT|AL_HUBER_LOSS);
alined.setLossScale(1.0);

// World Line endpoints (Don't need to be exact)
Eigen::Matrix<double,4,Eigen::Dynamic> X_w;

// Projected Line endpoints in normalized camera space
EIgen::Matrix<double,3,Eigen::Dynamic> x_c;

// Line-Based DLT
Eigen::Matrix4d pose = alined.poseFromLines(x_c,X_w);

// Iterative Approach
Eigen::Matrix4d tf = alined.poseFromLinesIterative(pose, x_c, X_w);
```


## Contributors

If you want to contribute to this library, feel free to send your contributor request to the mail address provided in the package.xml.

## License

 ALineD is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ALineD is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ALineD. If not, see <http://www.gnu.org/licenses/>.

