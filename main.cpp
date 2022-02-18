//essas bibliotecas estão em anexo e devem estar presentes no projeto final
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "GLutils.h"
#include "svpng.inc"
#include <iostream>

//usei para testar durante o desenvolvimento do software
#define ALTURA 700
#define LARGURA 700
#define NumeroIteracoes 3


//essas variaveis devem ser declaradas pois são usadas na renderização das imagens
VAO vao;
GLBuffer vbo,ebo;
ShaderProgram shaderProgram;

//usado para fazer interpolação entre dois pontos
vec2 Interpolacao(vec2 p1, vec2 p2, float ind){
    vec2 res=ind*p2+(1-ind)*p1;
    return res;
}

//função que acha o terceiro ponto de quadrado equilátero
vec2 TerceiroPontoTriangulo(vec2 p1,vec2 p3){
    double s60 = sin(300 * M_PI / 180.0);
    double c60 = cos(300 * M_PI / 180.0);
    return {
            c60 * (p1[0] - p3[0]) - s60 * (p1[1] - p3[1]) + p3[0],
            s60 * (p1[0] - p3[0]) + c60 * (p1[1] - p3[1]) + p3[1]
        };
}

//função responsável pela geração dos pontos da curva de koch para renderizar usando linhas
void geraPontosCurvaKoch(std::vector<vec2> &Pontos){
   for(int i=0;i<Pontos.size();i=i+4){
        std::vector<vec2> temporario;
        vec2 p1,p3;
        if((i+1)!=Pontos.size()){
            p1=Interpolacao(Pontos[i],Pontos[i+1],0.333333333333333333);
            p3=Interpolacao(Pontos[i],Pontos[i+1],0.666666666666666666);
        }else{
            p1=Interpolacao(Pontos[i],Pontos[0],0.333333333333333333);
            p3=Interpolacao(Pontos[i],Pontos[0],0.666666666666666666);
        }
        vec2 p2 = TerceiroPontoTriangulo(p1,p3);
        temporario.push_back(p1);
        temporario.push_back(p2);
        temporario.push_back(p3);
        auto iter = Pontos.insert(Pontos.begin() + (i+1), temporario.begin(), temporario.end());
   }
}

//função que gera os quadrados do tapete de sierpinski recursivamente
void GeraPontosSierpinski(vec2 ponto1,vec2 ponto2,int iteracaoAtual,int numeroIteracao,std::vector<vec2> &Pontos,std::vector<unsigned int> &indices){
    if(iteracaoAtual==numeroIteracao){
        int dif=(ponto2[0]-ponto1[0])/3,dif2=(ponto2[1]-ponto1[1])/3;
        Pontos.push_back({ponto1[0]+dif,ponto1[1]+dif2});
        indices.push_back(indices.size());
        Pontos.push_back({ponto1[0]+dif,ponto1[1]+dif2*2});
        indices.push_back(indices.size());
        Pontos.push_back({ponto1[0]+dif*2,ponto1[1]+dif2*2});
        indices.push_back(indices.size());
        Pontos.push_back({ponto1[0]+dif*2,ponto1[1]+dif2});
        indices.push_back(indices.size());
    }else{
        int dif=(ponto2[0]-ponto1[0])/3,dif2=(ponto2[1]-ponto1[1])/3;
        Pontos.push_back({ponto1[0]+dif,ponto1[1]+dif2});
        indices.push_back(indices.size());
        Pontos.push_back({ponto1[0]+dif,ponto1[1]+dif2*2});
        indices.push_back(indices.size());
        Pontos.push_back({ponto1[0]+dif*2,ponto1[1]+dif2*2});
        indices.push_back(indices.size());
        Pontos.push_back({ponto1[0]+dif*2,ponto1[1]+dif2});
        indices.push_back(indices.size());
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                if(i==1 && j==1)
                    continue;
                double temp=(ponto2[0]-ponto1[0])/3,temp2=(ponto2[1]-ponto1[1])/3;
                GeraPontosSierpinski({ponto1[0]+temp*i,ponto1[1]+temp2*j},{ponto1[0]+temp*(i+1),ponto1[1]+temp2*(j+1)},
                                     iteracaoAtual+1,numeroIteracao,Pontos,indices);
            }
        }
    }
}

//função principal que é responsável pela geração da curva de koch
void CurvaKoch(int largura,int altura,int numeroIteracao,std::vector<std::vector<int>> &MatrizCores){
	glewInit();

	shaderProgram = ShaderProgram{
		Shader{"SimpleShader.vert", GL_VERTEX_SHADER},
		Shader{"SimpleShader.frag", GL_FRAGMENT_SHADER}
	};
	glUseProgram(shaderProgram);

	vec2 p1={largura*0.2,altura*0.25},p3={largura*0.8,altura*0.25},p2=TerceiroPontoTriangulo(p1,p3);

	std::vector<vec2> P = {  p1,  p2, p3 };

    std::vector<unsigned int> indices;

    for(int i=0;i<numeroIteracao;i++)
        geraPontosCurvaKoch(P);

    for(int j=0;j<P.size();j++)
        indices.push_back(j);

	vao = VAO{true};
	glBindVertexArray(vao);

	vbo = GLBuffer{GL_ARRAY_BUFFER};
	vbo.data(P, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);

	ebo = GLBuffer{GL_ELEMENT_ARRAY_BUFFER};
	ebo.data(indices, GL_STATIC_DRAW);

	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);
	Uniform{"M"} = orthogonal(0, w, 0, h, -1, 1);
	Uniform{"C"} = vec3{0,0,0};

	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBindVertexArray(vao);
	glDrawElements(GL_LINE_LOOP, 12*pow(4,NumeroIteracoes-1) , GL_UNSIGNED_INT, 0);


	unsigned char pixels[largura*altura*4];
	//unsigned char vetorReverso[largura*altura*4];
	unsigned r,g,b,a;
    glReadPixels( 0, 0,largura, altura, GL_RGBA, GL_UNSIGNED_BYTE, pixels );

/*    int k=0;


    for(int i=tam-1; i>=0; i=i-4){
        a=pixels[i];
        b=pixels[i-1];
        g=pixels[i-2];
        r=pixels[i-3];;
        vetorReverso[k++] = r;
        vetorReverso[k++] = g;
        vetorReverso[k++] = b;
        vetorReverso[k++] = a;
    }

    for(int j=0;j<tam-1;j++)
        pixels[j]=vetorReverso[j];

*/
    //aqui a imagem com a curva de koch é salva no mesmo diretorio do projeto
    FILE* fp = fopen("koch.png", "wb");
    svpng(fp, largura, altura, pixels, 1);
    fclose(fp);

    //aqui a matriz de cores é preeenchida
    for(int i=0;i<largura;i++){
        std::vector<int> temp;
        for(int j=0;j<altura;j++){
            int index=4*(i * largura + j);
            if((pixels[index]
                +pixels[index+1]
                +pixels[index+2]
                +pixels[index+3])==(255*4)){
                //std::cout<<"0 ";
                temp.push_back(0);
            }else{
                //std::cout<<"1 ";
                temp.push_back(1);
            }
        }
        MatrizCores.push_back(temp);
        //std::cout<<'\n';
    }

    printf("\n Curva de koch gerada");
}

//função principal que gera o tapete de sierpinski
void TapeteSierpinski(int largura,int altura,int numeroIteracao,std::vector<std::vector<int>> &MatrizCores){
	glewInit();

	shaderProgram = ShaderProgram{
		Shader{"SimpleShader.vert", GL_VERTEX_SHADER},
		Shader{"SimpleShader.frag", GL_FRAGMENT_SHADER}
	};
	glUseProgram(shaderProgram);

	std::vector<vec2> P ,PTotal;

    std::vector<unsigned int> indicesTotal;

    for(int h=0;h<numeroIteracao;h++)
        GeraPontosSierpinski({0,0},{largura,altura},1,3,PTotal,indicesTotal);

    glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	for(int ind=0;ind<indicesTotal.size();){
        P.clear();
        P.push_back(PTotal[ind++]);
        P.push_back(PTotal[ind++]);
        P.push_back(PTotal[ind++]);
        P.push_back(PTotal[ind++]);

        vao = VAO{true};
        glBindVertexArray(vao);

        vbo = GLBuffer{GL_ARRAY_BUFFER};
        vbo.data(P, GL_STATIC_DRAW);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);

        ebo = GLBuffer{GL_ELEMENT_ARRAY_BUFFER};
        ebo.data(indicesTotal, GL_STATIC_DRAW);

        int w = glutGet(GLUT_WINDOW_WIDTH);
        int h = glutGet(GLUT_WINDOW_HEIGHT);
        Uniform{"M"} = orthogonal(0, w, 0, h, -1, 1);
        Uniform{"C"} = vec3{0,0,0};

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao);
        glDrawElements(GL_TRIANGLE_FAN, 4 , GL_UNSIGNED_INT, 0);

	}

	unsigned char pixels[largura*altura*4];
	//unsigned char vetorReverso[largura*altura*4];
	unsigned r,g,b,a;
    glReadPixels( 0, 0,largura, altura, GL_RGBA, GL_UNSIGNED_BYTE, pixels );

/*    int k=0;


    for(int i=tam-1; i>=0; i=i-4){
        a=pixels[i];
        b=pixels[i-1];
        g=pixels[i-2];
        r=pixels[i-3];;
        vetorReverso[k++] = r;
        vetorReverso[k++] = g;
        vetorReverso[k++] = b;
        vetorReverso[k++] = a;
    }

    for(int j=0;j<tam-1;j++)
        pixels[j]=vetorReverso[j];

*/
    //aqui a imagem com a curva de koch é salva no mesmo diretorio do projeto
    FILE* fp = fopen("Sierpinski.png", "wb");
    svpng(fp, largura, altura, pixels, 1);
    fclose(fp);

    //aqui a matriz de cores é preenchida
    for(int i=0;i<largura;i++){
        std::vector<int> temp;
        for(int j=0;j<altura;j++){
            int index=4*(i * largura + j);
            if((pixels[index]
                +pixels[index+1]
                +pixels[index+2]
                +pixels[index+3])==(255*4)){
                //std::cout<<"0 ";
                temp.push_back(0);
            }else{
                //std::cout<<"1 ";
                temp.push_back(1);
            }
        }
        MatrizCores.push_back(temp);
        //std::cout<<'\n';
    }

    printf("\n Tapete de Sierpinski gerado");

}

//main usados apenas como exemplo para a utilização das funções
int main(int argc, char* argv[]){
    //esse bloco de funções devem ser inseridos em algum lugar antes das funções que geram as imagens
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_MULTISAMPLE);
	glutInitWindowSize(LARGURA,ALTURA);
	glutInitContextVersion(3, 3);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutCreateWindow("Gerador de imagens");

	printf("GL Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL Version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

	//abaixo eu crio uma matriz responsável por armazenar as cores da imagem e chamo a função que gera a imagem desejada
	std::vector<std::vector<int>> MatrizCores;
	CurvaKoch(LARGURA,ALTURA,NumeroIteracoes,MatrizCores);
	//a parte abaixo pode ser usada para visualizar a matriz de cores
	/*for(int i=0;i<LARGURA;i++){
        for(int j=0;j<ALTURA;j++){
            std::cout<<MatrizCores[i][j]<<' ';
        }
        std::cout<<'\n';
	}*/
    MatrizCores.clear();

	TapeteSierpinski(LARGURA,ALTURA,NumeroIteracoes,MatrizCores);
	/*for(int i=0;i<LARGURA;i++){
        for(int j=0;j<ALTURA;j++){
            std::cout<<MatrizCores[i][j]<<' ';
        }
        std::cout<<'\n';
	}*/
}
