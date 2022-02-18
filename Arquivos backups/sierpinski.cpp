#include <GL/glew.h>
#include <GL/freeglut.h>
#include "GLutils.h"
#include "svpng.inc"
#include <iostream>

#define ALTURA 40
#define LARGURA 40
#define NumeroIteracoes 1

VAO vao;
GLBuffer vbo,ebo;
ShaderProgram shaderProgram;

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

int main(int argc, char **argv){
    glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_MULTISAMPLE);
	glutInitWindowSize(LARGURA,ALTURA);
	glutInitContextVersion(3, 3);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutCreateWindow("janela");
	printf("GL Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL Version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
	//init(LARGURA,ALTURA,NumeroIteracoes);


	std::vector<std::vector<int>> MatrizCores;
	TapeteSierpinski(LARGURA,ALTURA,NumeroIteracoes,MatrizCores);
	for(int i=0;i<LARGURA;i++){
        for(int j=0;j<ALTURA;j++){
            std::cout<<MatrizCores[i][j]<<' ';
        }
        std::cout<<'\n';
	}
    //glutMainLoop();
}
