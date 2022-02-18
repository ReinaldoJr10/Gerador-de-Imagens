#ifndef BEZIER_H
#define BEZIER_H

#include "vec.h"

inline vec2 cubic_bezier(const vec2 P[4], float t){
	vec2 Q[3] = {
		(1-t)*P[0] + t*P[1],
		(1-t)*P[1] + t*P[2],
		(1-t)*P[2] + t*P[3]
	};

	vec2 R[2] = {
		(1-t)*Q[0] + t*Q[1],
		(1-t)*Q[1] + t*Q[2],
	};

	return (1-t)*R[0] + t*R[1];
}

inline std::vector<vec2> sample_cubic_bezier(const vec2 P[4], int N){
	std::vector<vec2> C(N);
	for(int i = 0; i < N; i++){
		float t = i/(N-1.0);
		C[i] = cubic_bezier(P, t);
	}
	return C;
}

inline std::vector<vec2> sample_cubic_bezier_spline(const std::vector<vec2>& P, int N){
	std::vector<vec2> C;
	int n = P.size();
	for(int i = 0; i+3 < n; i+=3){
		std::vector<vec2> Q = sample_cubic_bezier(P.data() + i, N);
		C.insert(C.end(), Q.begin(), Q.end());
	}
	return C;
}

inline vec2 bezier(const std::vector<vec2>& P, float t){
	int n = P.size();
	if(n == 1)
		return P[0];

	std::vector<vec2> Q(n-1);
	for(int i = 0; i < n-1; i++)
		Q[i] = (1-t)*P[i] + t*P[i+1];

	return bezier(Q, t);
}

inline std::vector<vec2> sample_bezier(const std::vector<vec2>& P, int N){
	std::vector<vec2> C(N);
	for(int i = 0; i < N; i++){
		float t = i/(N-1.0);
		C[i] = bezier(P, t);
	}
	return C;
}

inline vec3 cubic_bezier(const vec3 P[4], float t){
	vec3 Q[3] = {
		(1-t)*P[0] + t*P[1],
		(1-t)*P[1] + t*P[2],
		(1-t)*P[2] + t*P[3]
	};

	vec3 R[2] = {
		(1-t)*Q[0] + t*Q[1],
		(1-t)*Q[1] + t*Q[2],
	};

	return (1-t)*R[0] + t*R[1];
}

inline std::vector<vec3> sample_cubic_bezier(const vec3 P[4], int N){
	std::vector<vec3> C(N);
	for(int i = 0; i < N; i++){
		float t = i/(N-1.0);
		C[i] = cubic_bezier(P, t);
	}
	return C;
}


#endif
