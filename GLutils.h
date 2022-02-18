#ifndef GLSL_UTILS_H
#define GLSL_UTILS_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <array>
#include <functional>

#ifndef COLOR_H
#define COLOR_H

#include <math.h>
#ifndef VEC_H
#define VEC_H

#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include <cmath>
#include <vector>

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


typedef unsigned char Byte;

struct Color{
	Byte r, g, b;
};

const Color white = {255, 255, 255};
const Color red = {255, 0, 0};
const Color blue = {0, 0, 255};
const Color green = {0, 255, 0};
const Color black = {0, 0, 0};
const Color cyan = {0, 255, 255};
const Color yellow = {255, 255, 0};
const Color magenta = {255, 0, 255};
const Color gray = {128, 128, 128};
const Color orange = {255, 165, 0};
const Color purple = {128, 0, 128};
const Color brown = {160, 128, 96};

inline float ByteToFloat(Byte v){
	return v/255.0f;
}

inline float clamp(float v, float a, float b){
	if(v < a)
		return a;
	else if(v > b)
		return b;
	return v;
}

inline Byte floatToByte(float v){
	return round(255*clamp(v, 0.0f, 1.0f));
}

inline vec3 toVec3(Color C){
	return {
		ByteToFloat(C.r),
		ByteToFloat(C.g),
		ByteToFloat(C.b)
	};
}

inline Color toColor(vec3 v){
	return {
		floatToByte(v(0)),
		floatToByte(v(1)),
		floatToByte(v(2))
	};
}

inline Color lerp(float t, Color A, Color B){
	vec3 vA = toVec3(A);
	vec3 vB = toVec3(B);
	return toColor(lerp(t, vA, vB));
}

inline Color bilerp(float u, float v, Color A, Color B, Color C, Color D){
	vec3 vA = toVec3(A);
	vec3 vB = toVec3(B);
	vec3 vC = toVec3(C);
	vec3 vD = toVec3(D);
	return toColor(bilerp(u, v, vA, vB, vC, vD));
}

#endif


////////////////////////////////////////////////////////////////////
void enable_debug();

////////////////////////////////////////////////////////////////////
inline int currentProgram(){
	int program;
	glGetIntegerv(GL_CURRENT_PROGRAM, &program);
	return program;
}

////////////////////////////////////////////////////////////////////

struct Uniform{
	int loc = -1;

	Uniform() = default;

	Uniform(std::string var){
		loc = glGetUniformLocation(currentProgram(), var.c_str());
	}

	Uniform& operator=(int v){
		glUniform1i(loc, v);
		return *this;
	}

	Uniform& operator=(float v){
		glUniform1f(loc, v);
		return *this;
	}

	Uniform& operator=(vec2 v){
		glUniform2fv(loc, 1, v.v);
		return *this;
	}

	Uniform& operator=(vec3 v){
		glUniform3fv(loc, 1, v.v);
		return *this;
	}

	Uniform& operator=(vec4 v){
		glUniform4fv(loc, 1, v.v);
		return *this;
	}

	Uniform& operator=(mat2 M){
		glUniformMatrix2fv(loc, 1, true, M.v);
		return *this;
	}

	Uniform& operator=(mat3 M){
		glUniformMatrix3fv(loc, 1, true, M.v);
		return *this;
	}

	Uniform& operator=(mat4 M){
		glUniformMatrix4fv(loc, 1, true, M.v);
		return *this;
	}

	Uniform& operator=(const std::vector<mat4>& M){
		glUniformMatrix4fv(loc, M.size(), true, (float*)M.data());
		return *this;
	}

	Uniform& operator=(const std::vector<float>& v){
		glUniform1fv(loc, v.size(), v.data());
		return *this;
	}

	template<size_t N>
	Uniform& operator=(const std::array<vec3, N>& v){
		glUniform3fv(loc, v.size(), v[0].v);
		return *this;
	}
};

////////////////////////////////////////////////////////////////////
struct UintResource{
	unsigned int id = 0;
	std::function<void(unsigned int)> deleter;

	operator unsigned int() const{ return id; }

	UintResource() = default;

	~UintResource(){
		if(id && deleter)
			deleter(id);
	}


	UintResource(UintResource&& other){
		id = other.id;
		other.id = 0;
	}

	UintResource& operator=(UintResource&& other){
		if(&other != this){
			UintResource tmp{std::move(other)};
			std::swap(id, tmp.id);
		}
		return *this;
	}
};

////////////////////////////////////////////////////////////////////
struct Shader : public UintResource{
	GLenum type;
	Shader() = default;
	Shader(std::string file, GLenum shader_type);
	void check();
};

////////////////////////////////////////////////////////////////////
struct ShaderProgram : public UintResource{
	ShaderProgram() = default;

	template<class...T>
	ShaderProgram(const T&... shaders){
		id = glCreateProgram();
		deleter = glDeleteProgram;
		init(shaders...);
	}

	int getAttribLocation(std::string attrib){
		return glGetAttribLocation(id, attrib.c_str());
	}
	private:
	void init(){
		glLinkProgram(id);
		check();
	}
	template<class...T>
	void init(const Shader& s, const T&...tail){
		glAttachShader(id, s);
		init(tail...);
		glDetachShader(id, s);
	}


	void check();
};

////////////////////////////////////////////////////////////////////
struct GLBuffer : public UintResource{
	GLenum type;

	GLBuffer() = default;

	GLBuffer(GLenum type) : type{type}{
		glGenBuffers(1, &id);
		deleter = [](unsigned int id){
			glDeleteBuffers(1, &id);
		};
	}

	template<class T>
	void data(const std::vector<T>& V, GLenum usage){
		glBindBuffer(type, id);
		glBufferData(type, V.size()*sizeof(T), V.data(), usage);
	}
};

////////////////////////////////////////////////////////////////////
struct VAO : public UintResource{
	VAO(bool init = false){
		if(!init)
			return;

		glGenVertexArrays(1, &id);

		deleter = [](unsigned int id){
			glDeleteVertexArrays(1, &id);
		};
	}
};

////////////////////////////////////////////////////////////////////
struct GLTexture : public UintResource{
	GLenum target;

	GLTexture() = default;

	GLTexture(GLenum target) : target{target}{
		glGenTextures(1, &id);
		glBindTexture(target, id);
		deleter = [](unsigned int id){
			glDeleteTextures(1, &id);
		};
	}

	void load(std::string filename, GLenum target = 0);
};

#endif

#ifndef TRANSFORMS_H
#define TRANSFORMS_H
#define _USE_MATH_DEFINES
#define M_PI           3.14159265358979323846

#include <cmath>

inline mat4 loadIdentity(){
    return Id<4>();
}

inline mat4 translate(vec3 u){
	return {
		1, 0, 0, u(0),
		0, 1, 0, u(1),
		0, 0, 1, u(2),
		0, 0, 0,    1
	};
}

inline mat4 translate(float x, float y, float z){
    return translate(vec3{x, y, z});
}

inline mat4 scale(float x, float y, float z){
	return {
		x, 0, 0, 0,
		0, y, 0, 0,
		0, 0, z, 0,
		0, 0, 0, 1
	};
}

inline mat4 rotate_x(float angle){
	float c = cos(angle);
	float s = sin(angle);

	return {
		1, 0,  0, 0,
		0, c, -s, 0,
		0, s,  c, 0,
		0, 0,  0, 1
	};
}

inline mat4 rotate_y(float angle){
	float c = cos(angle);
	float s = sin(angle);

	return {
		 c, 0, s, 0,
		 0, 1, 0, 0,
		-s, 0, c, 0,
		 0, 0, 0, 1
	};
}

inline mat4 rotate_z(float angle){
	float c = cos(angle);
	float s = sin(angle);

	return {
		c, -s, 0, 0,
		s,  c, 0, 0,
		0,  0, 1, 0,
		0,  0, 0, 1
	};
}

inline vec3 rodrigues(vec3 eixo,vec3 uni,float theta){
  //vec3 uni={1,1,1};
  uni=normalize(uni);
  return (dot( eixo,uni))*(1-cos(theta))*uni+
                  cos(theta)*eixo+
                  sin(theta)*cross(uni,eixo);
}

inline mat4 rotate(vec3 n, float theta){
  vec3 eixo1={1,0,0};
  vec3 eixo2={0,1,0};
  vec3 eixo3={0,0,1};

  vec3 rodrigues1= rodrigues(eixo1,n,theta);
  vec3 rodrigues2= rodrigues(eixo2,n,theta);
  vec3 rodrigues3= rodrigues(eixo3,n,theta);

	return {
    rodrigues1[0],rodrigues2[0],rodrigues3[0],0,
    rodrigues1[1],rodrigues2[1],rodrigues3[1],0,
    rodrigues1[2],rodrigues2[2],rodrigues3[2],0,
    0,0,0,1
  };
}

inline mat4 lookAt(vec3 eye, vec3 center, vec3 up){
	vec3 upNormalizado=normalize(up);
  vec3 f=normalize((eye-center));
  f=-f;
  vec3 s=normalize(cross(f,upNormalizado));
  vec3 u=cross(s,f);

  mat4 t1={
    s[0],s[1],s[2],0,
    u[0],u[1],u[2],0,
    -f[0],-f[1],-f[2],0,
    0,0,0,1
  };
  mat4 t2={
    1,0,0,-eye[0],
    0,1,0,-eye[1],
    0,0,1,-eye[2],
    0,0,0,1
  };

	return t1*t2;
}

inline mat4 orthogonal(float l, float r, float b, float t, float n, float f){
	return {
		2/(r-l),      0,     0,      (l+r)/(l-r),
			0,  2/(t-b),     0,      (b+t)/(b-t),
			0,        0, -2/(f-n),   (n+f)/(n-f),
			0,        0,      0,               1
	};
}

inline mat4 frustum(float l, float r, float b, float t, float n, float f){
	mat4 Persp = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0,-(n+f)/n, -f,
		0, 0, -1/n, 0
	};
	mat4 Volume = {
		2/(r-l),      0,     0,      (l+r)/(l-r),
			0,  2/(t-b),     0,      (b+t)/(b-t),
			0,        0, 2/(f-n),    (n+f)/(n-f),
			0,        0,      0,               1
	};
	return Volume*Persp;
}


inline mat4 perspective(float fovy, float aspect, float Near, float Far){
    float t,b,r,l;
    t=Near*tan((fovy*M_PI/180)/2);
    b=-1*t;
    r=t*aspect;
    l=-r;
    return frustum(l,r,b,t,Near,Far);
}

#endif
