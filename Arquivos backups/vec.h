#ifndef VEC_H
#define VEC_H

#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include <array>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>

template<int M, int N=M>
struct Mat{
    float v[M*N];

    float& operator()(int i, int j=0){
        return v[i*N + j];
    }

    float operator()(int i, int j=0) const{
        return v[i*N + j];
    }

    Mat operator+(const Mat& B) const{
        Mat R;
        for(int i = 0; i < M*N; i++)
            R.v[i] = v[i] + B.v[i];
        return R;
    }

    Mat operator-(const Mat& B) const{
        Mat R;
        for(int i = 0; i < M*N; i++)
            R.v[i] = v[i] - B.v[i];
        return R;
    }

    Mat operator-() const{
        Mat R;
        for(int i = 0; i < M*N; i++)
            R.v[i] = -v[i];
        return R;
    }

    template<int P>
    Mat<M, P> operator*(const Mat<N, P>& B) const{
        const Mat& A = *this;
        Mat<M, P> R;
        for(int i = 0; i < M; i++)
            for(int j = 0; j < P; j++){
                R(i, j) = 0;
                for(int k = 0; k < N; k++)
                    R(i, j) += A(i, k)*B(k, j);
            }
        return R;
    }

    Mat<M,1> get_column(int j) const{
        const Mat& A = *this;
        Mat<M,1> C;
        for(int i = 0; i < M; i++)
            C(i) = A(i, j);
        return C;
    }

    float& operator[](int j){
        return v[j];
    }

    float operator[](int j)const{
        return v[j];
    }

    void set_column(const Mat<M,1>& C, int j){
        Mat& A = *this;

        for(int i = 0; i < M; i++)
            A(i, j) = C(i);
    }

    friend void print(const Mat& A){
        for(int i = 0; i < M; i++){
            for(int j = 0; j < N; j++)
                printf("%10.4f ", A(i, j));
            printf("\n");
        }
    }
    friend void print(std::string name, const Mat& A){
        printf("%s = ", name.c_str());
        if(M > 1)
            puts("");
        print(A);
    }

    friend Mat<N, M> transpose(const Mat& A){
        Mat<N, M> r;
        for(int i = 0; i < M; i++)
            for(int j = 0; j < N; j++)
                r(j, i) = A(i, j);
        return r;
    }

    friend Mat operator*(float x, const Mat& A){
        Mat r = A;
        for(float& v: r.v)
            v *= x;
        return r;
    }

    friend float norm2(const Mat& A){
        float r = 0;
        for(float x: A.v)
            r += x*x;
        return r;
    }

    friend float norm(const Mat& A){
        return sqrt(norm2(A));
    }
};

using mat2 = Mat<2>;
using mat3 = Mat<3>;
using mat4 = Mat<4>;

using mat2x2 = Mat<2,2>;
using mat2x3 = Mat<2,3>;
using mat2x4 = Mat<2,4>;

using mat3x2 = Mat<3,2>;
using mat3x3 = Mat<3,3>;
using mat3x4 = Mat<3,4>;

using mat4x2 = Mat<4,2>;
using mat4x3 = Mat<4,3>;
using mat4x4 = Mat<4,4>;

template<int M>
Mat<M> Id(){
    Mat<M> R;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++)
            R(i, j) = (i==j);
    return R;
}

template<int M>
Mat<M-1> submatrix(const Mat<M>& A, int i0, int j0){
    Mat<M-1> R;
    int k = 0;
    for(int i = 0; i < M; i++){
        if(i == i0)
            continue;

        for(int j = 0; j < M; j++){
            if(j == j0)
                continue;

            R.v[k] = A(i, j);
            k++;
        }
    }

    return R;
}

template<int M>
float determinant(const Mat<M>& A);

template<int M>
float cofactor(const Mat<M>& A, int i, int j){
    float r = determinant(submatrix(A, i, j));
    if((i+j)%2!=0)
       r = -r;
    return r;
}

inline float determinant(const Mat<1>& A){
    return A(0, 0);
}

inline float determinant(const mat2& A){
    return A(0, 0)*A(1,1) - A(1, 0)*A(0, 1);
}

inline float determinant(const mat3& A){
    return A(0, 0)*(A(1,1)*A(2,2) - A(1, 2)*A(2,1))
          -A(1, 0)*(A(0,1)*A(2,2) - A(0, 2)*A(2, 1))
          +A(2, 0)*(A(0,1)*A(1,2) - A(0, 2)*A(1,1));
}

template<int M>
float determinant(const Mat<M>& A){
    float det = 0;
    for(int i = 0; i < M; i++)
        det += A(i, 0)*cofactor(A, i, 0);
    return det;
}

template<int M>
Mat<M> cofactor_matrix(const Mat<M>& A){
    Mat<M> R;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++)
            R(i,j) = cofactor(A, i, j);
    return R;
}

template<int M>
Mat<M> inverse(const Mat<M>& A){
    float det = determinant(A);
    if(fabs(det) < 1e-5)
        return {NAN};

    Mat<M> adjoint = transpose(cofactor_matrix(A));

    return 1/det*adjoint;
}

inline mat4 toMat4(mat3 A){
	return {
        A(0, 0), A(0, 1), A(0, 2), 0,
        A(1, 0), A(1, 1), A(1, 2), 0,
        A(2, 0), A(2, 1), A(2, 2), 0,
		      0,       0,       0, 1
	};
}

inline mat3 toMat3(mat4 A){
	return {
        A(0, 0), A(0, 1), A(0, 2),
        A(1, 0), A(1, 1), A(1, 2),
        A(2, 0), A(2, 1), A(2, 2)
	};
}

#endif

template<int M>
using Vec = Mat<M, 1>;

using vec2 = Vec<2>;
using vec3 = Vec<3>;
using vec4 = Vec<4>;

template<int M>
Vec<M> operator*(const Vec<M>& U, const Vec<M>& V){
    Vec<M> R;
    for(int i = 0; i < M; i++)
        R(i) = U(i)*V(i);
    return R;
}

template<int M>
Vec<M> operator/(const Vec<M>& U, const Vec<M>& V){
    Vec<M> R;
    for(int i = 0; i < M; i++)
        R(i) = U(i)/V(i);
    return R;
}

template<int M>
float dot(const Vec<M>& U, const Vec<M>& V){
    float r = 0;
    for(int i = 0; i < M; i++)
        r += U(i)*V(i);
    return r;
}

template<int M>
Vec<M> normalize(const Vec<M>& u){
	return (1.0/norm(u))*u;
}

template<int M, int N>
std::vector<Vec<M>> operator*(const Mat<M, N>& A, const std::vector<Vec<N>>& P){
	int n = P.size();
	std::vector<Vec<M>> Q(n);
	for(int i = 0; i < n; i++)
		Q[i] = A*P[i];
	return Q;
}

inline vec2 toVec2(vec2 u){
    return u;
}

inline vec2 toVec2(vec3 u){
    return {u(0), u(1)};
}

inline vec2 toVec2(vec4 u){
    return {u(0), u(1)};
}

inline vec3 toVec3(vec2 u, float z=0){
    return {u(0), u(1), z};
}

inline vec3 toVec3(vec3 u){
    return u;
}

inline vec3 toVec3(vec4 u){
    return {u(0), u(1), u(2)};
}

inline vec4 toVec4(vec2 u, float z=0, float w=1){
    return {u(0), u(1), z, w};
}

inline vec4 toVec4(vec3 u, float w=1){
    return {u(0), u(1), u(2), w};
}

inline vec4 toVec4(vec4 u){
    return u;
}

inline std::vector<vec2> loadCurve(std::string filename){
	std::ifstream in(filename);
	int n = 0;
	in >> n;
	std::vector<vec2> P(n);
	for(vec2& v: P)
		in >> v(0) >> v(1);
	return P;
}

template<int M>
Vec<M> lerp(float t, Vec<M> A, Vec<M> B){
	return (1-t)*A + t*B;
}

template<int M>
Vec<M> bilerp(float s, float t, Vec<M> A, Vec<M> B, Vec<M> C, Vec<M> D){
	Vec<M> E = lerp(s, A, B);
	Vec<M> F = lerp(s, C, D);
	return lerp(t, E, F);
}

inline float find_mix_param(vec2 v, vec2 v0, vec2 v1){
	vec2 d = v1 - v0;
	return dot(d, v-v0)/dot(d,d);
}

inline vec3 cross(vec3 u, vec3 v){
	return {
		 u(1)*v(2) - u(2)*v(1),
		 u(2)*v(0) - u(0)*v(2),
		 u(0)*v(1) - u(1)*v(0)
	};
}

template<int M, class...V, int N = 1+sizeof...(V)>
Mat<M, N> toMat(const Vec<M>& C, V... tail){
	std::array<Vec<M>, N> columns = {C, tail...};
	Mat<M, N> R;
	for(int j = 0; j < N; j++)
		R.set_column(columns[j], j);
    return R;
}

/*****************************************************************************/
/* TRANSFORMAÇÕES EM COORDENADAS HOMOGÊNEAS */
inline vec2 operator*(const mat3& A, vec2 P){
	vec3 Q = A*toVec3(P, 1);
	return 1/Q(2)*toVec2(Q);
}

inline std::vector<vec2> operator*(const mat3& A, const std::vector<vec2>& Ps){
	std::vector<vec2> Q;

	for(vec2 P: Ps)
		Q.push_back( A*P );

	return Q;
}
/*****************************************************************************/


#endif
