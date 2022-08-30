library(shinythemes)
library(plotly)
library(dplyr)
library(shinyalert)
library(shinyWidgets)
library(DT)
library(ggrepel)
library(tibble)
library(ggplot2)
library(shiny)
library(shinydashboard)
library(agricolae)
tabelat=data.frame(1:50,
                   qt(1-0.25,1:50),
                   qt(1-0.10,1:50),
                   qt(1-0.05,1:50),
                   qt(1-0.025,1:50),
                   qt(1-0.01,1:50),
                   qt(1-0.005,1:50),
                   qt(1-0.0025,1:50),
                   qt(1-0.001,1:50),
                   qt(1-0.0005,1:50))
colnames(tabelat)=c("GL",0.25,0.10,0.05,0.025,0.01,0.005,0.0025,0.001,0.0005)
tabelaz=round(t(matrix(pnorm(seq(0,3.99,0.01))-0.5,ncol=40)),4)
tabelaz=data.frame(cbind(seq(0,3.9,0.1),tabelaz))
colnames(tabelaz)=c("Z0",0:9)
tabelaf=function(sig=0.05){
  var2=c(1:30,40,50,60,80,100,120,240)
  tabf=list()
  for(i in 1:length(var2))
    tabf[[i]]=qf(1-sig,c(1:15,20,40,60,120,240),df2 = var2[i])
  tabf=data.frame(round(t(matrix(unlist(tabf),ncol=37)),3))
  tabf=cbind(var2,tabf)
  colnames(tabf)=c(" ",1:15,20,40,60,120,240)
  tabf
}

plot_normal = function(start, end, p, z, syb="Z",inf=NA,sup=NA){
  mean.1 <- 0
  sd.1 <- 1
  zstart <- 2
  zend <- 3
  zcritical <- 1.65
  my_col <- "#00998a"
  x <- seq(from = mean.1 - 3 * sd.1,
           to = mean.1 + 3 * sd.1,
           by = .01)
  MyDF <- data.frame(x = x, y = dnorm(x, mean = mean.1, sd = sd.1))
  shade_curve <-
    function(MyDF,
             zstart,
             zend,
             fill = "white",
             alpha = 1) {
      geom_area(
        data = subset(MyDF, x >= mean.1 + zstart * sd.1
                      & x < mean.1 + zend * sd.1),
        aes(y = y),
        fill = fill,
        color = NA,
        alpha = alpha
      )
    }
  p1a <- ggplot(MyDF, aes(x = x, y = y)) +
    geom_line(color="blue",size=1)+
    geom_polygon(fill = "red", color = "blue",alpha=0.5) +
    shade_curve(
      MyDF = MyDF,
      zstart = start,
      zend = end,
      fill = "white"
    ) + 
    labs(title=parse(text = paste(syb,"[cal] == ",z,"~~p-value == ", p,sep = "")))+
    theme_minimal() + 
    theme(axis.text = element_text(size=12))+
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(lty=2,xintercept = c(inf,sup))+
    annotate(geom = "segment",x=0,xend = 0,y=0,yend=max(MyDF$y),lty=2,size=0.7)+
    ylab("") +
    xlab("")
  p1a
}
testezp=function(x1,
                 x2,
                 n1,
                 n2,
                 alpha=0.05,
                 alternative="greater"){
  f1=x1/n1
  f2=x2/n2
  p=(x1+x2)/(n1+n2)
  zcal=(f1-f2)/sqrt((p*(1-p))*(1/n1+1/n2))
  if(alternative=="greater"){
    ztab=qnorm(1-alpha/2)
    zpvalue=(1-pnorm(abs(zcal),lower.tail = TRUE))*2
    p1=plot_normal(-abs(zcal),abs(zcal),zpvalue,zcal,syb = "Z",inf = -abs(ztab),sup = abs(ztab))}
  if(alternative=="esquerda"){
    ztab=qnorm(alpha)
    zpvalue=pnorm(zcal,lower.tail = TRUE)
    p1=plot_normal(zcal,Inf,zpvalue,zcal,syb = "Z",inf = ztab)}
  if(alternative=="direita"){
    ztab=qnorm(1-alpha)
    zpvalue=(1-pnorm(zcal,lower.tail = TRUE))
    p1=plot_normal(-Inf,zcal,zpvalue,zcal,syb = "Z",sup=ztab)}
  list(f1=f1,
       f2=f2,
       p=p,
       "Statistic" = zcal,
       "Z_standard" = ztab,
       "normal_plot" = p1,
       "p-value" = zpvalue)
}
testez=function(x1,
                x2,
                var1,
                var2,
                n1,
                n2,
                alpha=0.05,
                alternative="greater"){
  zcal=(x1-x2)/sqrt(var1^2/n1+var2^2/n2)
  if(alternative=="greater"){
    ztab=qnorm(1-alpha/2)
    zpvalue=(1-pnorm(abs(zcal),lower.tail = TRUE))*2
    p1=plot_normal(-abs(zcal),abs(zcal),zpvalue,zcal,syb = "Z",inf = -abs(ztab),sup = abs(ztab))}
  if(alternative=="esquerda"){
    ztab=qnorm(alpha)
    zpvalue=pnorm(zcal,lower.tail = TRUE)
    p1=plot_normal(zcal,Inf,zpvalue,zcal,syb = "Z",inf = ztab)}
  if(alternative=="direita"){
    ztab=qnorm(1-alpha)
    zpvalue=(1-pnorm(zcal,lower.tail = TRUE))
    p1=plot_normal(-Inf,zcal,zpvalue,zcal,syb = "Z",sup = ztab)}
  list("Statistic" = zcal,
       "Z_standard" = ztab,
       "normal_plot" = p1,
       "p-value" = zpvalue)
}
testet_vari=function(x1,
                     x2,
                     s1,
                     s2,
                     n1,
                     n2,
                     alpha=0.05,
                     alternative="greater"){
  sc=((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)
  tcal=(x1-x2)/(sqrt(sc)*sqrt(1/n1+1/n2))
  if(alternative=="greater"){
    ttab=qt(1-alpha/2,df = n1+n2-2)
    tpvalue=(1-pt(abs(tcal),df =n1+n2-2,lower.tail = TRUE))*2
    p1=plot_normal(-abs(tcal),abs(tcal),tpvalue,tcal,syb = "t",inf = -abs(ttab),sup = abs(ttab))}
  if(alternative=="esquerda"){
    ttab=qt(alpha,df = n1+n2-2)
    tpvalue=pt(tcal,df =n1+n2-2,lower.tail = TRUE)
    p1=plot_normal(tcal,Inf,tpvalue,tcal,syb = "t",inf = ttab)}
  if(alternative=="direita"){
    ttab=qt(1-alpha,df = n1+n2-2)
    tpvalue=(1-pt(tcal,df =n1+n2-2, lower.tail = TRUE))
    p1=plot_normal(-Inf,tcal,tpvalue,tcal,syb = "t",sup = ttab)}
  list("Statistic" = tcal,
       "t_value" = ttab,
       "t_plot" = p1,
       "var_c"=sc,
       "p-value" = tpvalue)
}

testet_vard=function(x1,
                     x2,
                     s1,
                     s2,
                     n1,
                     n2,
                     alpha=0.05,
                     alternative="greater"){
  tcal=(x1-x2)/sqrt(s1^2/n1+s2^2/n2)
  v=(s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n2-1)+(s2^2/n2)^2/(n2-1)) 
  if(alternative=="greater"){
    ttab=qt(1-alpha/2,df = v)
    tpvalue=(1-pt(abs(tcal),df = v,lower.tail = TRUE))*2
    p1=plot_normal(-abs(tcal),abs(tcal),tpvalue,tcal,syb = "t",inf = -abs(ttab),sup = abs(ttab))}
  if(alternative=="esquerda"){
    ttab=qt(alpha,df = v)
    tpvalue=pt(tcal,df = v,lower.tail = TRUE)
    p1=plot_normal(tcal,Inf,tpvalue,tcal,syb = "t",inf = ttab)}
  if(alternative=="direita"){
    ttab=qt(1-alpha,df = n1+n2-2)
    tpvalue=(1-pt(tcal,df = v, lower.tail = TRUE))
    p1=plot_normal(-Inf,tcal,tpvalue,tcal,syb = "t",sup = ttab)}
  list("Statistic" = tcal,
       "v"=v,
       "t_value" = ttab,
       "t_plot" = p1,
       "p-value" = tpvalue)
}


testefhomog=function(s1,s2,n1,n2,sig=0.05){
  s=c(s1,s2)
  n=c(n1,n2)
  f=max(s)/min(s)
  gl1=n[s==max(s)]-1
  gl2=n[s==min(s)]-1
  ft1=qf(1-sig/2,gl1,gl2)
  ft2=qf(sig/2,gl1,gl2)
  pvalor=(1-pf(f,gl1,gl2))*2
  list(nume=max(s),
       den=min(s),
       "fmax"=f,
       "ft1"=ft1,
       "ft2"=ft2,
       gl1=gl1,
       gl2=gl2,
       p=pvalor)
}

testez_one=function(x, mu, var, n, alpha=0.05, alternative="greater"){
  zcal=(x-mu)/(var/sqrt(n))
  if(alternative=="greater"){
    ztab=qnorm(1-alpha/2)
    zpvalue=(1-pnorm(abs(zcal),lower.tail = TRUE))*2
    p1=plot_normal(-abs(zcal),abs(zcal),zpvalue,zcal,syb = "Z",inf = -abs(ztab),sup = abs(ztab))}
  if(alternative=="esquerda"){
    ztab=qnorm(alpha)
    zpvalue=pnorm(zcal,lower.tail = TRUE)
    p1=plot_normal(zcal,Inf,zpvalue,zcal,syb = "Z",inf = ztab)}
  if(alternative=="direita"){
    ztab=qnorm(1-alpha)
    zpvalue=(1-pnorm(zcal,lower.tail = TRUE))
    p1=plot_normal(-Inf,zcal,zpvalue,zcal,syb = "Z",sup = ztab)}
  list("Statistic" = zcal,
       "Z_standard" = ztab,
       "normal_plot" = p1,
       "p-value" = zpvalue)}

testet_one=function(x, mu, var, n, alpha=0.05, alternative="greater"){
  tcal=(x-mu)/(var/sqrt(n))
  if(alternative=="greater"){
    ttab=qt(1-alpha/2,df = n-1)
    tpvalue=(1-pt(abs(tcal),df =n-1,lower.tail = TRUE))*2
    p1=plot_normal(-abs(tcal),abs(tcal),tpvalue,tcal,syb = "t",inf = -abs(ttab),sup = abs(ttab))}
  if(alternative=="esquerda"){
    ttab=qt(alpha,df = n-1)
    tpvalue=pt(tcal,df =n-1,lower.tail = TRUE)
    p1=plot_normal(tcal,Inf,tpvalue,tcal,syb = "t",inf = ttab)}
  if(alternative=="direita"){
    ttab=qt(1-alpha,df = n-1)
    tpvalue=(1-pt(tcal,df =n-1, lower.tail = TRUE))
    p1=plot_normal(-Inf,tcal,tpvalue,tcal,syb = "t",sup = ttab)}
  list("Statistic" = tcal,
       "t_value" = ttab,
       "t_plot" = p1,
       "p-value" = tpvalue)}

####################################################################
ui <- dashboardPage(
  dashboardHeader(title='R for education'),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home",tabName = "home",icon = icon("home")),br(),
      menuItem("  Análise exploratória de dados",
      menuSubItem("Análise descritiva",tabName = "descritiva"),
      menuSubItem("Média Ponderada",tabName = "mediapond"),
      menuSubItem("Gráfico de caixas",tabName = "boxplot")),
      menuItem(" Tabela de frequência",
               menuSubItem("Qualitativo",tabName = "tabquali"),
               menuSubItem("Quantitativo discreto",tabName="tabquantidis"),
               menuSubItem("Quantitativo contínuo",tabName="tabquanticont")),
      menuItem("  Teste de hipóteses",
               menuSubItem('Teste Z para uma média', tabName='taba'),
      menuSubItem('Teste t para uma média', tabName='tabb'),
      menuSubItem('Teste Z para uma proporção', tabName='tabc'),
      menuSubItem('Teste Z para duas médias', tabName='tab1a'),
      menuSubItem('Teste Z para duas proporções', tabName='tab1b'),
      menuItem('Teste t para duas médias', tabName='tabs',
               menuSubItem('Var. desc. \ne igual', tabName='tab2a'),
               menuSubItem('Var. desc. \ne desigual', tabName='tab2b')),
      menuSubItem('Teste F para duas variâncias', tabName='tab1var'),
      menuSubItem('Teste t pareado', tabName='tab2c')),
      menuItem("Experimentação agrícola",
      menuItem('Normalidade dos erros', tabName='norm'),
      menuItem('Homogeneidade das variâncias', tabName='homog'),
      menuItem("Anova",tabName = "anova",
               menuSubItem('DIC', tabName='anova1'),
               menuSubItem('DBC', tabName='anova2'),
               menuSubItem('DQL', tabName='anova3'),
               menuSubItem('DIC em esquema fatorial', tabName='anova4'),
               menuSubItem('DBC em esquema fatorial', tabName='anova5'))),
      
      menuItem("Correlação e regressão",
      menuItem("Correlação Pearson",tabName = "corre"),
      menuItem("Regressão linear simples",tabName = "regre")),
      
      menuItem("Tabelas",
      menuItem("Tabela t",tabName = "tabelat"),
      menuItem("Tabela Z",tabName = "tabelaz"),
      menuItem("Tabela F",tabName = "tabelaf"))
    )
  ),
  dashboardBody(
    fluidRow(
      withMathJax(),
      tabItems(
        tabItem(tabName = "home",
                box(h2("Seja vem-vindo à nossa plataforma de ensino de estatística básica"),
                    img(src="https://github.com/AgronomiaR/images/blob/main/r_education.jpg?raw=true", width="100%",heigth="100%"),
                    width = 12),
                box(title = "Objetivo",
                    "Essa plataforma visa atender aos alunos de graduação, pós-graduação e todos os profissionais que usam e estudam estatística. Ela foi criado com o objetivo de auxiliar na compreensão de cálculos e conceitos básicos de estatística, permitindo ao usuário visualizar cada etapa do processo de análise. Essa plataforma foi construída inteiramente em linguagem R, usando o pacote shiny e shinydashboard.",width = 8),
                box(title = "Organização e desenvolvimento",
                    "Desenvolvedor: Prof. Msc. Gabriel Danilo Shimizu",br(),br(),
                    "Universidade Estadual de Londrina, Departamento de estatística",width = 4),br(),br()),
        tabItem(tabName = "descritiva",
                box(textInput("dadosdescritiva",
                              "Digitar dados (Separar números por vírgula e decimal por ponto)",
                              value = 0),
                    actionButton("run","Executar",width = "32%"),width = 12),
                box(title = "Média aritmética", 
                    "A média é definida como a somatória das observações dividida pelo número de observações, assim temos que:",br(),
                    withMathJax(helpText("$$\\bar{x} = \\frac{x_1 + x_2 + ... + x_3 + x_4}{n}$$")),
                    uiOutput("media"),width = 12),
                box(title = "Mediana", 
                    "É o valor que divide o conjunto ordenado de valores em duas partes com igual número de elementos, ou seja, 50% das observações ficam acima da mediana e 50% ficam abaixo.",br(),
                    uiOutput("mediana"),
                    textOutput("mediana1"),
                    uiOutput("mediana2"),
                    width = 12),
                box(title = "Variância amostral",
                    "Variância amostral é  a soma dos quadrados dos desvios, dividida pelo total de observações menos um.",br(),
                    withMathJax(helpText("$$S^2 = \\frac{\\sum(x_i-\\bar{x})^2}{n-1}$$")),
                    uiOutput("variancia"),width = 12),
                box(title = "Desvio-padrão amostral",
                    "O desvio-padrão amostral é a raiz quadrada da variância e recebe a mesma unidade de medida da variável", br(),
                    withMathJax(helpText("$$S = \\sqrt{S^2}$$")),
                    uiOutput("desviopadrao"),width = 12),
                box(title = "Coeficiente de variação",
                    "O coeficiente de variação é a razão entre o desvio-padrão e a média. É uma medida relativa que independe da unidade de medida",br(),
                    withMathJax(helpText("$$CV (\\%) = \\frac{S}{\\bar{x}}\\times100$$")),
                    uiOutput("cv"),width = 12)),
        
        
        tabItem(tabName = "tabquali",
                box(textInput("dadosqualitab",
                              "Digitar dados (Separar por vírgula)",
                              value = 0),
                    actionButton("gotabquali",div(icon("play-circle"),"  RUN"),width = "32%"),width = 12),
                box(width = 12,
                    shiny::dataTableOutput("tabqualifreq")),
                box(width = 6,
                    plotlyOutput("plotpizza")),
                box(width = 6,
                    plotlyOutput("plotqualifreq"))),
        
        tabItem(tabName = "tabquantidis",
                h1("EM CONSTRUÇÃO"),
                box(textInput("dadosquantidis",
                              "Digitar dados (Separar por vírgula)",
                              value = 0),
                    actionButton("gotabquantdis",div(icon("play-circle"),"  RUN"),width = "32%"),width = 12),
                box(width = 12,
                    shiny::dataTableOutput("tabquantidisfreq")),
                box(width = 6,
                    plotlyOutput("plotquantidis"))),
        
        tabItem(tabName = "tabquanticont",
                h1("EM CONSTRUÇÃO"),
                box(textInput("dadosquanticont",
                              "Digitar dados (Separar por vírgula)",
                              value = 0),
                    actionButton("gotabquanticont",div(icon("play-circle"),"  RUN"),width = "32%"),width = 12),
                box(width = 12,
                    "Número de classes:",br(),
                    withMathJax(helpText("$$Método \\ 1: \\sqrt{n}$$")),
                    uiOutput("k1"),
                    withMathJax(helpText("$$Método \\ 2: \\sqrt{n} - 1$$")),
                    uiOutput("k2"),
                    withMathJax(helpText("$$Método \\ 3: 1 + 3,3 \\times log(n)$$")),
                    uiOutput("k3")),
                box(width = 12,
                    numericInput("k_definir",label = "Número de classes (k)",value = 5),
                    "Amplitude:",
                    uiOutput("amplitude"),br(),
                    "Intervalo de classe calculado:",
                    uiOutput("intervalo_classe"),br(),
                    actionButton("gotabquanticont1",div(icon("play-circle"),"  RUN"),width = "32%"),br(),br(),
                    "Definir intervalo de classe:",
                    numericInput("int_definir",label = "Intervalo de classes (h)",value = 0),br(),
                    h3("Gerar Tabela:"),
                    actionButton("gotabquanticont2",div(icon("play-circle"),"  RUN"),width = "32%"),br(),br(),
                    shiny::dataTableOutput("tabquantcontfreq")),
                box(width = 12,
                    plotlyOutput("plotquantfreq"))),

        tabItem(tabName='mediapond',
                h2("Média ponderada"),
                box(fileInput('uploadedcsvpond', "File",accept = c('.csv','.xls','.xlsx'),
                              width = "200%",multiple = TRUE),
                    "Obs. Importar uma planilha em excel ou csv. Deve conter apenas a primeira coluna com os pesos e a segunda com as observações. A primeira linha como cabeçalho"),
                box(numericInput("planilhapond","Aba",value = 1,width = "50%"),
                    selectInput("deccsvpond","Separação CSV",choices = c(",",";","."),selected = ",", width = "50%"),
                    actionButton("gopond",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Média Ponderada",
                    "A média ponderada é a somatória das observações multiplicadas pelos respectivos pesos e dividida pela somatória dos pesos, conforme a seguir:", br(),
                    withMathJax(helpText("$$\\bar{x} = \\frac{\\sum^N_{i=1}\\bar{x}_i n_i}{n}$$")),
                    uiOutput("mediaponderada"),width = 12)),
        tabItem(tabName='boxplot',
                h2("Gráfico de caixas"),
                box(textInput("dadoscaixas",
                              "Digitar dados (Separar números por vírgula e decimal por ponto)",
                              value = 0),
                    actionButton("gobox",div(icon("play-circle"),"  RUN"),width = "32%"),width = 12),
                box(title = "Ordenar conjunto de dados",
                    "A primeira etapa é ordenar o conjunto de dados de forma crescente",br(),
                    uiOutput("box1")),
                box(title = "Mediana",
                    "A primeira medida de posição a ser encontrada é a mediana, que dividi em 50% o conjunto de dados",br(),
                    uiOutput("box2")),
                box(title = "Primeiro quartil (Q1)",
                    "A segunda medida de posição é o primeiro quartil que representa os primeiros 25% do conjunto de dados. Para encontrar, a forma mais fácil é calcular a mediana dos valores abaixo da mediana.",br(),
                    uiOutput("box3")),
                box(title = "Terceiro quartil (Q3)",
                    "A terceira medida de posição é o terceiro quartil que representa os primeiros 75% do conjunto de dados. Para encontrar, a forma mais fácil é calcular a mediana dos valores acima da mediana.",br(),
                    uiOutput("box4")),
                box(title = "Limite superior",
                    "O limite superior é calculado pela fórmula:",br(),
                    withMathJax(helpText("$$LS = Q3 + 1,5(Q3-Q1)$$")),
                    "Se o valor encontrado foi maior que o maior valor do conjunto de dados, o limite superior é substituido pelo máximo",br(),
                    uiOutput("box5"),
                    textOutput("box5a")),
                box(title = "Limite inferior",
                    "O limite inferior é calculado pela fórmula:",br(),
                    withMathJax(helpText("$$LI = Q1 - 1,5(Q3-Q1)$$")),
                    "Se o valor encontrado foi menor que o menor valor do conjunto de dados, o limite inferior é substituido pelo mínimo",br(),
                    uiOutput("box6"),
                    textOutput("box6a")),
                box(title = "Outliers",
                    "Se houver observações abaixo ou acima do limite inferior e superior, respectivamente, é considerado um outlier.",br(),
                    br(),"Conclui-se que:",br(),
                    textOutput("box7")),
                box(title = "Gráfico",
                    plotlyOutput("plot_caixas"))),
        tabItem(tabName='taba',
                h2("Teste Z para uma média populacional com variância conhecida"),
                box(title ="Introdução", 
                    "Teste Z para uma média populacional (μ) com desvios-padrão populacional conhecido (σ).",br(),br(),
                    "Selecione as hipóteses nula e alternativa, digite o nível de significância, a média da amostra, o desvio-padrão da população, e o tamanho amostral", br(),br(),
                    selectInput("m1h0z1","Hipótese Nula (H0)",choices = c("=","≥","≤"),selected = "="),
                    selectInput("m1h1z1","Hipótese Alternativa (H1)",choices = c("≠","<",">"),selected = "≠"),br(),br(),
                    uiOutput("m1hipz1"),
                    uiOutput("m1hipz2")),
                box(title="Gráfico",
                    plotOutput(outputId = "m1testez2")),
                box(title="Inserção de dados",
                    box(numericInput("m1z21","Média Amostral",value = 10),
                        numericInput("m1zp21","Média Populacional",value = 10),
                        numericInput("m1sdz21","Desvio-padrão Populacional",value = 2)),
                    box(numericInput("m1nz21","Tamanho amostral",value = 30)),
                    box(numericInput("m1sigz21","Nível de significância",value = 0.05)),
                    actionButton("m1goz2",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    withMathJax(helpText(
                      paste("$$Z=\\frac{\\bar{x} - \\mu}{\\frac{\\sigma}{\\sqrt{n}}}$$"))),
                    uiOutput("m1contavalorzz"),br(),br(),
                    box(title = "Valor de p",
                        textOutput("m1conclpz2")),
                    box(title = "Valor de Z tabelado",
                        textOutput("m1ztabe2")))),
        tabItem(tabName='tabb',
                h2("Teste t para uma média populacional com variância desconhecida"),
                box(title ="Introdução", 
                    "Teste t para uma média populacional (μ) com desvios-padrão populacional desconhecido (σ).",br(),br(),
                    "Selecione as hipóteses nula e alternativa, digite o nível de significância, a média da amostra, o desvio-padrão amostral, e o tamanho amostral", br(),br(),
                    selectInput("m2h0z1","Hipótese Nula (H0)",choices = c("=","≥","≤"),selected = "="),
                    selectInput("m2h1z1","Hipótese Alternativa (H1)",choices = c("≠","<",">"),selected = "≠"),br(),br(),
                    uiOutput("m2hipz1"),
                    uiOutput("m2hipz2")),
                box(title="Gráfico",
                    plotOutput(outputId = "m2testez2")),
                box(title="Inserção de dados",
                    box(numericInput("m2z21","Média Amostral",value = 10),
                        numericInput("m2zp21","Média Populacional",value = 10),
                        numericInput("m2sdz21","Desvio-padrão amostral",value = 2)),
                    box(numericInput("m2nz21","Tamanho amostral",value = 30)),
                    box(numericInput("m2sigz21","Nível de significância",value = 0.05)),
                    actionButton("m2goz2",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    withMathJax(helpText(
                      paste("$$t=\\frac{\\bar{x} - \\mu}{\\frac{S}{\\sqrt{n}}}$$"))),
                    uiOutput("m2contavalorzz"),
                    box(title = "Valor de p",
                        textOutput("m2conclpz2")),
                    box(title = "Valor de Z tabelado",
                        textOutput("m2ztabe2")))),
        tabItem(tabName='tabc',
                h2("Teste Z para uma proporção populacional"),
                box("Em construção !!! ")),
        tabItem(tabName='tab1a',
                h2( "Teste Z para duas médias populacionais com variância conhecida"),
                box(title ="Introdução", 
                    "Teste Z para duas médias populacionais (μ1 e μ2) com desvios-padrões populacionais conhecidos (σ1 e σ2).",br(),br(),
                    "Selecione as hipóteses nula e alternativa, digite o nível de significância, as médias da amostra, os desvios-padrão da população, 
              os tamanhos das amostras", br(),br(),
                    selectInput("h0z2","Hipótese Nula (H0)",choices = c("=","≥","≤"),selected = "="),
                    selectInput("h1z2","Hipótese Alternativa (H1)",choices = c("≠","<",">"),selected = "≠"),br(),br(),
                    uiOutput("hip1"),
                    uiOutput("hip2")),
                box(title="Inserção de dados",
                    box(numericInput("mediasz21","Medias Amostral 1",value = 10),
                        numericInput("mediasz22","Medias Amostral 2",value = 10)),
                    box(numericInput("sdz21","Desvio-padrão Populacional 1",value = 2),
                        numericInput("sdz22","Desvio-padrão Populacional 2",value = 2)),
                    box(numericInput("nz21","Tamanho amostral 1",value = 30),
                        numericInput("nz22","Tamanho amostral 2",value = 30)),
                    box(numericInput("sigz21","Nível de significância",value = 0.05)),
                    actionButton("goz2",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    withMathJax(helpText(
                      paste("$$Z=\\frac{\\bar{x_1} - \\bar{x_2}}{\\sqrt{\\frac{\\sigma_1^2}{n_1}+\\frac{\\sigma_2^2}{n_2}}}$$"))),
                    uiOutput("contavalorzz"),
                    box(title = "Valor de p",
                        textOutput("conclpz2")),
                    box(title = "Valor de Z tabelado",
                        textOutput("ztabe2"))),
                box(title="Gráfico",
                    plotOutput(outputId = "testez2"))
        ),
        
        tabItem(tabName='tab1b',
                h2( "Teste Z para duas proporções populacionais com variância conhecida"),
                box(title ="Introdução", 
                    "Teste Z para duas proporções populacionais (p1 e p2).",br(),br(),
                    "Selecione as hipóteses nula e alternativa, digite o nível de significância, o número de ocorrências e os tamanhos das amostras", br(),br(),
                    selectInput("h0zp2","Hipótese Nula (H0)",choices = c("=","≥","≤"),selected = "="),
                    selectInput("h1zp2","Hipótese Alternativa (H1)",choices = c("≠","<",">"),selected = "≠"),
                    uiOutput("hipp1"),
                    uiOutput("hipp2")),
                box(title="Inserção de dados",
                    box(numericInput("pz21","Total ocorrências 1",value = 10),
                        numericInput("pz22","Total ocorrências 2",value = 10)),
                    box(numericInput("pnz21","Tamanho amostral 1",value = 30),
                        numericInput("pnz22","Tamanho amostral 2",value = 30)),
                    box(numericInput("psigz21","Nível de significância",value = 0.05)),
                    actionButton("pgoz2",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    withMathJax(helpText(
                      paste("$$f_1=\\frac{x_1}{n_1}$$"))),
                    uiOutput("f1p"),
                    withMathJax(helpText(
                      paste("$$f_2=\\frac{x_2}{n_2}$$"))),
                    uiOutput("f2p"),
                    withMathJax(helpText(
                      paste("$$p=\\frac{x_1+x_2}{n_1+x_2}$$"))),
                    uiOutput("pp"),
                    withMathJax(helpText(
                      paste("$$Z=\\frac{f_1 - f_2}{p(1-p)\\sqrt{\\frac{1}{n_1}+\\frac{1}{n_2}}}$$"))),
                    uiOutput("pcontaz"),
                    box(title = "Valor de p",
                        textOutput("pconcl2")),
                    box(title = "Valor de Z tabelado",
                        textOutput("pztabe2"))),
                box(title="Gráfico",
                    plotOutput(outputId = "ptestez2"))
        ),
        
        tabItem(tabName='tab2a',
                h2( "Teste t para duas médias populacionais com variância desconhecida e igual"),
                box(title ="Introdução", 
                    "Teste t para duas médias populacionais (μ1 e μ2) com desvios-padrões desconhecidos mas iguais.",br(),br(),
                    "Selecione as hipóteses nula e alternativa, digite o nível de significância, as médias da amostra, os desvios-padrão da amostra e 
              os tamanhos das amostras", br(),br(),
                    selectInput("h0t2","Hipótese Nula (H0)",choices = c("=","≥","≤"),selected = "="),
                    selectInput("h1t2","Hipótese Alternativa (H1)",choices = c("≠","<",">"),selected = "≠"),br(),br(),
                    uiOutput("hipt1"),
                    uiOutput("hipt2")),
                box(title="Inserção de dados",
                    box(numericInput("mediast21","Medias Amostral 1",value = 10),
                        numericInput("mediast22","Medias Amostral 2",value = 10)),
                    box(numericInput("sdt21","Desvio-padrão amostral 1",value = 2),
                        numericInput("sdt22","Desvio-padrão amostral 2",value = 2)),
                    box(numericInput("nt21","Tamanho amostral 1",value = 30),
                        numericInput("nt22","Tamanho amostral 2",value = 30)),
                    box(numericInput("sigt21","Nível de significância",value = 0.05)),
                    actionButton("got2",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    "Variância conjunta:",br(),
                    withMathJax(helpText(
                      paste("$$S_c^2=\\frac{(n_1-1)S_1^2+(n_2-1)S_2^2}{n_1+n_2-2}$$"))),
                    uiOutput("scvalortt"),
                    "Estatística do teste:",br(),
                    withMathJax(helpText(
                      paste("$$t=\\frac{\\bar{x_1} - \\bar{x_2}}{S_c\\sqrt{\\frac{1}{n_1}+\\frac{1}{n_2}}}$$"))),
                    uiOutput("contavalortt"),
                    
                    box(title = "Valor de p",
                        textOutput("conclpt2")),
                    box(title = "Valor de t tabelado",
                        textOutput("ttabe2"))),
                box(title="Gráfico",
                    plotOutput(outputId = "testetg2"))
        ),
        
        tabItem(tabName='tab2b',
                h2( "Teste t para duas médias populacionais com variância desconhecida e desigual"),
                box(title ="Introdução", 
                    "Teste t para duas médias populacionais (μ1 e μ2) com desvios-padrões desconhecidos e desiguais.",br(),br(),
                    "Selecione as hipóteses nula e alternativa, digite o nível de significância, as médias da amostra, os desvios-padrão da amostra e 
              os tamanhos das amostras", br(),br(),
                    selectInput("h0tt2","Hipótese Nula (H0)",choices = c("=","≥","≤"),selected = "="),
                    selectInput("h1tt2","Hipótese Alternativa (H1)",choices = c("≠","<",">"),selected = "≠"),br(),br(),
                    uiOutput("hiptt1"),
                    uiOutput("hiptt2")),
                box(title="Inserção de dados",
                    box(numericInput("mediastt21","Medias Amostral 1",value = 10),
                        numericInput("mediastt22","Medias Amostral 2",value = 10)),
                    box(numericInput("sdtt21","Desvio-padrão amostral 1",value = 2),
                        numericInput("sdtt22","Desvio-padrão amostral 2",value = 2)),
                    box(numericInput("ntt21","Tamanho amostral 1",value = 30),
                        numericInput("ntt22","Tamanho amostral 2",value = 30)),
                    box(numericInput("sigtt21","Nível de significância",value = 0.05)),
                    actionButton("gott2",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    "Estatística do Teste (t):",br(),
                    withMathJax(helpText(
                      paste("$$Z=\\frac{\\bar{x_1} - \\bar{x_2}}{\\sqrt{\\frac{S_1^2}{n_1}+\\frac{S_2^2}{n_2}}}$$"))),
                    uiOutput("contavalorttt"),
                    "Graus de liberdade de Welch-Satterthwaite (v):",br(),
                    withMathJax(helpText(
                      paste("$$v = \\frac{(\\frac{S_1^2}{n_1}+\\frac{S_2^2}{n_2})^2}{\\frac{(S_1^2/n_1)^2}{n_1-1}+\\frac{(S_2^2/n_2)^2}{n_2-1}}$$"))),
                    uiOutput("vvalorttt"),
                    box(title = "Valor de p",
                        textOutput("conclptt2")),
                    box(title = "Valor de t tabelado",
                        textOutput("tttabe2"))),
                box(title="Gráfico",
                    plotOutput(outputId = "testettg2"))
        ),
        tabItem(tabName='tab1var',
                h2("Teste F para duas variâncias populacionais"),
                box(title ="Introdução",
                    "Teste F para igualdade de duas variâncias populacionais",
                    "Digite o nível de significância, as variâncias das amostras e o tamanho amostral", br(),br(),
                    withMathJax(helpText("$$H_0: \\ Variâncias \\ são \\ iguais$$")),
                    withMathJax(helpText("$$H_1: \\ Variâncias \\ são \\ diferentes$$"))),
                box(box(numericInput("var21","Variância amostral 1",value = 1),
                        numericInput("var22","Variância amostral 2",value = 2)),
                    box(numericInput("nvar21","Tamanho amostral 1",value = 30),
                        numericInput("nvar22","Tamanho amostral 2",value = 30)),
                    box(numericInput("sigvar","Nível de significância",value = 0.05)),
                    actionButton("govar2",div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Estatística do teste e valor de p",
                    "Estatística do Teste (F):",br(),
                    uiOutput("contavalorvar"),
                    "Graus de liberdade 1:",br(),
                    uiOutput("GL1"),
                    "Graus de liberdade 2:",br(),
                    uiOutput("GL2"),
                    box(title = "Valor de p",
                        textOutput("conclpvar")),
                    box(title = "Valor de F tabelado",
                        textOutput("vartabe")))),
        tabItem(tabName='tab2c',
                h2("Teste t pareado para duas amostras dependentes"),
                box(textInput("pareadoamostra1","Digitar dados Amostra 1 (Separar números por vírgula e decimal por ponto)",value = 0),
                    textInput("pareadoamostra2","Digitar dados Amostra 2 (Separar números por vírgula e decimal por ponto)",value = 0),
                    actionButton("gobox",div(icon("play-circle"),"  RUN"),width = "32%"),width = 12),
                box("Em construção !!! ")),
        tabItem(tabName='anova1',
                h2("Delineamento inteiramente casualizado"),
                "Obs. Use apenas para dados balanceados, ou seja, com mesmo número de repetições por tratamento",br(),
                box(fileInput('uploadedcsvaov', "File",accept = c('.csv','.xls','.xlsx'),width = "200%",multiple = TRUE),
                    numericInput("planilhaaov","Aba",value = 1,width = "50%"),
                    selectInput("deccsvaov","Separação CSV",choices = c(",",";","."),selected = ",", width = "50%"),
                    actionButton("goaov",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Entrada de dados",
                    box(selectInput("tratA", "Tratamentos", choices = "Tratamento" ),width = "50%",background = "navy"),
                    
                    box(selectInput("respA", "Variável resposta", choices = "Resposta" ),width = "50%",background = "navy")),br(),
                box("Correção (C)",
                    uiOutput("correcao1"),
                    "Graus de liberdade (GL)",
                    uiOutput("gltotal"),
                    uiOutput("gltrat"),
                    uiOutput("glresiduo"),
                    "Soma de quadrados total (SQT)",
                    uiOutput("sqt"),
                    "Soma de quadrados tratamento (SQTrat)",
                    uiOutput("sqtrat"),
                    "Soma de quadrados do resíduo",
                    uiOutput("sqres"),
                    "Quadrado médio do tratamento",
                    uiOutput("qmtrat"),
                    "Quadrado médio do resíduo",
                    uiOutput("qmres"),
                    "F calculado",
                    uiOutput("fcal"),width = 12),
                box("Resumo do quadro da Anova",
                    shiny::dataTableOutput("quadroanovadic"),width = 12),
                box("Teste de Tukey",
                    numericInput("alphatukeydic","Nível de significância",value = 0.05),
                    uiOutput("dmstukeydic"),
                    shiny::dataTableOutput("tukeytabela"),
                    width = 12)),
        tabItem(tabName='anova2',
                h2("Delineamento em blocos casualizados completos"),
                "Obs. Use apenas para dados balanceados, ou seja, com mesmo número de repetições por tratamento",br(),
                box(fileInput('uploadedcsvaov1', "File",accept = c('.csv','.xls','.xlsx'),width = "200%",multiple = TRUE),
                    numericInput("planilhaaov1","Aba",value = 1,width = "50%"),
                    selectInput("deccsvaov1","Separação CSV",choices = c(",",";","."),selected = ",", width = "50%"),
                    actionButton("goaov1",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Entrada de dados",
                    box(selectInput("tratB", "Tratamentos", choices = "Tratamento" ),width = "50%",background = "navy"),
                    box(selectInput("blocoB", "Blocos", choices = "Bloco" ),width = "50%",background = "navy"),
                    box(selectInput("respB", "Variável resposta", choices = "Resposta" ),width = "50%",background = "navy")),br(),
                box("Correção (C)",
                    uiOutput("correcao2"),
                    "Graus de liberdade (GL)",
                    uiOutput("gltotal1"),
                    uiOutput("gltrat1"),
                    uiOutput("glbloc1"),
                    uiOutput("glresiduo1"),
                    "Soma de quadrados total (SQT)",
                    uiOutput("sqt1"),
                    "Soma de quadrados tratamento (SQTrat)",
                    uiOutput("sqtrat1"),
                    "Soma de quadrados blocos (SQB)",
                    uiOutput("sqbloc1"),
                    "Soma de quadrados do resíduo",
                    uiOutput("sqres1"),
                    "Quadrado médio do tratamento",
                    uiOutput("qmtrat1"),
                    "Quadrado médio do bloco",
                    uiOutput("qmbloc1"),
                    "Quadrado médio do resíduo",
                    uiOutput("qmres1"),
                    "F calculado",
                    uiOutput("fcal1"),
                    uiOutput("fcalbloco1"),width = 12),
                box("Resumo do quadro da Anova",
                    shiny::dataTableOutput("quadroanovadbc"),width = 12),
                box("Teste de Tukey",
                    numericInput("alphatukeydbc","Nível de significância",value = 0.05),
                    uiOutput("dmstukeydbc"),
                    shiny::dataTableOutput("tukeytabeladbc"),
                    width = 12)),
        tabItem(tabName='anova3',
                h2("Delineamento em quadrado latino"),
                box("Em construção !!! ")),
        tabItem(tabName='corre',
                h2("Análise de correlação linear de pearson"),
                box(fileInput('uploadedcsvcor', "File",accept = c('.csv','.xls','.xlsx'),
                              width = "200%",multiple = TRUE),
                    numericInput("planilhacor","Aba",value = 1,width = "50%"),
                    selectInput("deccsvcor","Separação CSV",choices = c(",",";","."),
                                selected = ",", width = "50%"),
                    actionButton("gocor",#"RUN",
                                 div(icon("play-circle"),"  RUN"), width = "100%")),
                box(title = "Entrada de dados",
                    box(selectInput("X", "Variável I", choices = "I" ),width = "50%",background = "navy"),
                    box(selectInput("Y", "Variável II", choices = "II" ),width = "50%",background = "navy")),br(),
                box("Média de X",br(),
                    uiOutput("mediaX"),
                    "Média de Y",br(),
                    uiOutput("mediaY"),
                    "Desvio-padrão de X",br(),
                    uiOutput("desvioX"),
                    "Desvio-padrão de Y",br(),
                    uiOutput("desvioY"),
                    "Coeficiente de correlação de Pearson",br(),
                    withMathJax(helpText("$$r = \\frac{\\sum^{n}_{i=1}(X_i-\\bar{X})(Y_i-\\bar{Y})}{(n-1)S_xS_y} $$")),
                    uiOutput("pearson"),
                    width = 12)),
        tabItem(tabName='regre',
                h2("Análise de regressão linear simples"),
                box("Em construção !!! ")),
        tabItem(tabName='tabelat',
                h2("Tabela t de Student unicaudal superior"),
                dataTableOutput("tabelatestet")),
        tabItem(tabName='tabelaz',
                h2("Tabela Z da Distribuição Normal Padrão"),
                dataTableOutput("tabelatestez")),
        tabItem(tabName='tabelaf',
                h2("Tabela F de Fischer-Snedecor"),br(),
                numericInput("sigF","Nível de significância",value = 0.05),
                "Coluna representa o grau de liberdade do numerador e linha representa grau de liberdade do denominador.",
                dataTableOutput("tabelatestef"))
      )
    )
  ))

server <- function(input, output, session) {
  
  act_data=reactiveVal(0)
    trocar=eventReactive(input$run,{
      as.vector(as.numeric(unlist(strsplit(input$dadosdescritiva,","))))
    })
    output$media=renderUI({
      a=paste(trocar(),c(rep("+",length(trocar())-1)," "),
              collapse = '')
      withMathJax(helpText(paste("$$\\bar{x} = \\frac{",a,"}{",length(trocar()),"} = ",
                                 mean(trocar()),"$$")))})
    output$mediana=renderUI({
      a=paste(sort(trocar()),
              c(rep("\\to",length(trocar())-1)," "),
              collapse = '')
      withMathJax(helpText(paste("$$Ordem \\ crescente = ",a,"$$")))})
    output$mediana1=renderText({
      is.even <- function(x) x %% 2 == 0
      teste=is.even(length(trocar()))
      if(teste==TRUE){"Como o número de observações é par, a mediana é a média dos dois valores centrais"}
      if(teste==FALSE){"Como o número de observações é ímpar, a mediana é exatamente o valor central"}
    })
    
    output$mediana2=renderUI({
      withMathJax(helpText(paste("$$Assim, \\ a \\ mediana \\ é \\ ",median(trocar()),"$$")))})
    output$variancia=renderUI({
      a=paste("(",trocar(),"-",
              round(mean(trocar()),3),")^2",
              c(rep("+",length(trocar())-1)," "),
              collapse = '')
      withMathJax(helpText(paste("$$S^2 = \\frac{",a,"}{",length(trocar()),"-1} = ",round(var(trocar()),3),"$$")))})
    output$desviopadrao=renderUI({
      a=paste("(",trocar(),"-",
              round(mean(trocar()),3),")^2",
              c(rep("+",length(trocar())),""),
              collapse = '')
      withMathJax(helpText(paste("$$S = \\sqrt{\\frac{",a,"}{",length(trocar()),"-1}} = ","\\sqrt{",round(var(trocar()),3),
                                 "} = ",round(sd(trocar()),3),"$$")))})
    output$cv=renderUI({
      withMathJax(helpText(paste("$$CV = \\frac{",round(sd(trocar()),3),"}{",round(mean(trocar()),3),"} \\times 100 = ",round(sd(trocar())/mean(trocar())*100,3) ,"$$")))})
  # })
  
  ###########################################################################################
  ## Media ponderada
  dados_pond <- eventReactive(input$gopond,{
    infile <- input$uploadedcsvpond
    if (is.null(infile)){return(cat("Selecione a coluna de peso e resposta"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhapond))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvpond)}
    peso=tbl[,1]
    resp=tbl[,2]
    list("resp"=resp,"peso"=peso)})
  output$mediaponderada = renderUI({
    resp=dados_pond()$resp
    peso=dados_pond()$peso
    a=paste(unlist(peso),"\\times",unlist(resp),c(rep("+",length(resp)-1),""),collapse = '')
    mediap=round(sum(resp*peso)/sum(peso),3)
    withMathJax(helpText(paste("$$\\bar{x} = \\frac{",a , "}{",sum(peso),"} = ",mediap,"$$")))})

  ###########################################################################################
  corre=eventReactive(input$gobox,{
    list("resp"=as.vector(as.numeric(unlist(strsplit(input$dadoscaixas,",")))))
  })
  
  output$plot_caixas=renderPlotly(plot_ly(y=corre()$resp,type = "box", quartilemethod="exclusive"))
  output$box1 = renderUI({
    a=paste(sort(corre()$resp),
            c(rep("\\to",length(corre()$resp)-1)," "),
            collapse = '')
    withMathJax(helpText(paste("$$Ordem \\ crescente = ",a,"$$")))})
  output$box2 = renderUI(withMathJax(helpText(paste("$$Mediana = ",median(corre()$resp),"$$"))))
  quart = function(x) {
    x <- sort(x)
    n <- length(x)
    m <- (n+1)/2
    if (floor(m) != m) {
      l <- m-1/2; u <- m+1/2
    } else {
      l <- m-1; u <- m+1
    }
    c(Q1=median(x[1:l]), Q3=median(x[u:n]))
  }
  output$box3 = renderUI({
    quartis=quart(corre()$resp)
    withMathJax(helpText(paste("$$Q1 = ", quartis[1],"$$")))})
  output$box4 = renderUI({
    quartis=quart(corre()$resp)
    withMathJax(helpText(paste("$$Q3 = ",quartis[2],"$$")))})
  output$box5 = renderUI({
    quartis=quart(corre()$resp)
    q1=quartis[1]
    q3=quartis[2]
    withMathJax(helpText(paste("$$LS = ",q3,"+ 1.5(",q3,"-",q1,") = ",q3+1.5*(q3-q1),"$$")))})
  output$box6 = renderUI({
    quartis=quart(corre()$resp)
    q1=quartis[1]
    q3=quartis[2]
    withMathJax(helpText(paste("$$LI = ",q1,"- 1.5(",q3,"-",q1,") = ",q1-1.5*(q3-q1),"$$")))})
  output$box5a = renderText({
    quartis=quart(corre()$resp)
    q1=quartis[1]
    q3=quartis[2]
    maximo=max(corre()$resp)
    ls=q3+1.5*(q3-q1)
    if(ls>maximo){paste("Como o limite superior é maior que o máximo, o LS é o maior valor do conjunto de dados. Assim LS = ",max(corre()$resp))}})
  output$box6a = renderText({
    quartis=quart(corre()$resp)
    q1=quartis[1]
    q3=quartis[2]
    minimo=max(corre()$resp)
    li=q1-1.5*(q3-q1)
    if(li<minimo){paste("Como o limite inferior é menor que o mínimo, o LI é o menor valor do conjunto de dados. Assim LI = ",min(corre()$resp))}
  })
  out <- eventReactive(input$gobox,{
    quartis=quart(corre()$resp)
    q1=quartis[1]
    q3=quartis[2]
    maximo=max(corre()$resp)
    minimo=max(corre()$resp)
    li=q1-1.5*(q3-q1)
    ls=q3+1.5*(q3-q1)
    if(ls>maximo){ls=maximo}
    if(li<minimo){li=minimo}
    nmax=corre()$resp[ls>corre()$resp]
    nmin=corre()$resp[li<corre()$resp]
    list("nmax"=nmax,"nmin"=nmin)})
  output$box7=renderText({
    if(length(out()$nmax)>0 | length(out()$nmin)>0){a = "Há pontos discrepantes"}
    if(length(out()$nmax)==0 | length(out()$nmin)==0){a = "Não há pontos discrepantes"}
    a
  })
  
  
  #========================================================================================================
  # tabela de frequencia variável qualitativa
  #========================================================================================================
  qualitativa=eventReactive(input$gotabquali,{
    list("resp"=as.vector(unlist(strsplit(input$dadosqualitab,","))))
  })
  output$tabqualifreq=shiny::renderDataTable({
    resp=qualitativa()$resp
    tabela=table(resp)
    tabela=data.frame(tabela)
    colnames(tabela)=c("Variável","Fa")
    tabela$Fr=tabela$Fa/sum(tabela$Fa)
    tabela$`Fr%`=tabela$Fa/sum(tabela$Fa)*100
    tabela$`Faa`=cumsum(tabela$Fa)
    tabela$`Fra`=cumsum(tabela$Fa/sum(tabela$Fa))
    tabela
    })
  output$plotqualifreq=renderPlotly({
    resp=qualitativa()$resp
    tabela=table(resp)
    tabela=data.frame(tabela)
    a=ggplot(tabela,aes(x=resp,y=Freq))+
      geom_col()+theme_bw()
    ggplotly(a)})
  output$plotpizza=renderPlotly({
    resp=qualitativa()$resp
    tabela=table(resp)
    tabela=data.frame(tabela)
    plot_ly(tabela,labels=tabela$resp,
            values=tabela$Freq,type="pie")})
  
  #========================================================================================================
  # tabela de frequencia variável quantitativo discreto
  #========================================================================================================
  quantitativodis=eventReactive(input$gotabquantdis,{
    list("resp"=as.vector(unlist(strsplit(input$dadosquantidis,","))))
  })
  output$tabquantidisfreq=shiny::renderDataTable({
    resp=quantitativodis()$resp
    tabela=table(resp)
    tabela=data.frame(tabela)
    colnames(tabela)=c("Variável","Fa")
    tabela$Fr=tabela$Fa/sum(tabela$Fa)
    tabela$`Fr%`=tabela$Fa/sum(tabela$Fa)*100
    tabela$`Faa`=cumsum(tabela$Fa)
    tabela$`Fra`=cumsum(tabela$Fa/sum(tabela$Fa))
    tabela
  })
  output$plotquantidis=renderPlotly({
    resp=quantitativodis()$resp
    tabela=table(resp)
    tabela=data.frame(tabela)
    a=ggplot(tabela,aes(x=resp,y=Freq))+
      geom_col()+theme_bw()
    ggplotly(a)})

  #========================================================================================================
  # tabela de frequencia variável quantitativo
  #========================================================================================================
  quanticont=eventReactive(input$gotabquanticont,{
    list("resp"=as.vector(as.numeric(unlist(strsplit(input$dadosquanticont,",")))))})
  output$k1=renderUI({withMathJax(helpText(paste("$$ k = \\sqrt{",length(quanticont()$resp),"} = ",sqrt(length(quanticont()$resp)),"$$")))})
  output$k2=renderUI({withMathJax(helpText(paste("$$ k = \\sqrt{",length(quanticont()$resp),"} - 1 = ",
                                                 sqrt(length(quanticont()$resp))-1,"$$")))})
  output$k3=renderUI({withMathJax(helpText(paste("$$ k = 1 + 3,3 log(",length(quanticont()$resp),") = ",1+3.3*log(length(quanticont()$resp)),"$$")))})
  
  quant_classes=eventReactive(input$gotabquanticont1,{
    resp=as.vector(as.numeric(unlist(strsplit(input$dadosquanticont,","))))
    amp=max(resp)-min(resp)
    k=input$k_definir
    list("amp"=amp,"k"=k,
         "resp"=resp)})
  output$amplitude=renderUI({
    withMathJax(helpText(paste("$$ Amplitude (Amp) = ",max(quant_classes()$resp),"-",min(quant_classes()$resp),
                               " = ",quant_classes()$amp,"$$")))})
  output$intervalo_classe=renderUI({
    withMathJax(helpText(paste("$$ Intervalo (h) = \\frac{",quant_classes()$amp,
                                                 "}{",quant_classes()$k,"} = ",
                               quant_classes()$amp/quant_classes()$k,"$$")))})
  quant_classes_final=eventReactive(input$gotabquanticont2,{
    resp=quanticont()$resp
    k=input$k_definir
    int=input$int_definir
    tabela=cut(sort(resp),seq(from=min(resp),by=int,length=k+1))
    data.frame(table(tabela))})
  output$tabquantcontfreq=shiny::renderDataTable({
    quant_classes_final()})
  plot_ly_quant=eventReactive(input$gotabquanticont2,{
    resp=quanticont()$resp
    k=input$k_definir
    int=input$int_definir
    tabela=cut(sort(resp),seq(min(resp),by=int,length=k+1))
    tabela=data.frame(table(tabela))
    tabela$x=seq(min(resp),by=int,length=k)+int/2
    a=ggplot(tabela,aes(x=x,y=Freq))+
      geom_col(width = int,fill="lightblue",color="black")+
      scale_x_continuous(breaks = seq(min(resp),by=int,length=k+1))
    ggplotly(a)})
  output$plotquantfreq=renderPlotly({plot_ly_quant()})
  
  #========================================================================================================
  # teste z para uma amostra
  #========================================================================================================
  saidaz1=eventReactive(input$m1goz2,{
    if(input$m1h1z1=="≠"){alternativa="greater"}
    if(input$m1h1z1=="<"){alternativa="esquerda"}
    if(input$m1h1z1==">"){alternativa="direita"}
    testez_one(x = as.numeric(input$m1z21),
               mu = as.numeric(input$m1zp21),
               var = as.numeric(input$m1sdz21),
               n = as.numeric(input$m1nz21),
               alpha = as.numeric(input$m1sigz21),
               alternative = alternativa)})
  output$m1testez2=renderPlot(saidaz1()$normal_plot)
  output$m1hip1=renderUI(withMathJax(helpText(paste("$$H_0: \\mu_1",input$m1h0z2,input$m1zp21,"$$"))))
  output$m1hip2=renderUI(withMathJax(helpText(paste("$$H_1: \\mu_1",input$m1h1z2,input$m1zp21,"$$"))))
  output$m1contavalorzz=renderUI({
    withMathJax(helpText(
      paste("$$Z_{cal}=\\frac{",input$m1z21,"-",input$m1zp21,"}{\\frac{",input$m1sdz21,"}{\\sqrt{",
            input$m1nz21,"}}",
            "} = ",saidaz1()$Statistic,"$$")))})
  output$m1conclpz2=renderText(paste("O valor p é de",saidaz1()$`p-value`, "e o nível de significância é de",input$m1sigz21,". Assim, conclue-se que:",
                                   ifelse(saidaz1()$`p-value`>input$m1sigz21,"Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$m1ztabe2 = renderText({paste("Com o nível de significância é ",input$m1sigz21, "o valor crítico para um teste de hipótese",
                                    ifelse(input$m1h1z1=="≠","bicaudal","unilateral"),"é Ztab = ",saidaz1()$Z_standard)
    paste("O limite crítico é",
          ifelse(input$m1h1z1=="≠",paste(-round(saidaz1()$Z_standard,4)," e ",
                                         round(saidaz1()$Z_standard,4)),
                 ifelse(input$m1h1z1==">",
                        paste(round(saidaz1()$Z_standard,4)),
                        paste(round(saidaz1()$Z_standard,4)))))})
  
  
  #========================================================================================================
  # teste t para uma amostra
  #========================================================================================================
  saidaz2=eventReactive(input$m2goz2,{
    if(input$m2h1z1=="≠"){alternativa="greater"}
    if(input$m2h1z1=="<"){alternativa="esquerda"}
    if(input$m2h1z1==">"){alternativa="direita"}
    testet_one(x = as.numeric(input$m2z21),
               mu = as.numeric(input$m2zp21),
               var = as.numeric(input$m2sdz21),
               n = as.numeric(input$m2nz21),
               alpha = as.numeric(input$m2sigz21),
               alternative = alternativa)})
  output$m2testez2=renderPlot(saidaz2()$t_plot)
  output$m2hip1=renderUI(withMathJax(helpText(paste("$$H_0: \\mu_1",input$m2h0z2,input$m2zp21,"$$"))))
  output$m2hip2=renderUI(withMathJax(helpText(paste("$$H_1: \\mu_1",input$m2h1z2,input$m2zp21,"$$"))))
  output$m2contavalorzz=renderUI({
    withMathJax(helpText(
      paste("$$t_{cal}=\\frac{",input$m2z21,"-",input$m2zp21,"}{\\frac{",input$m2sdz21,"}{\\sqrt{",
            input$m2nz21,"}}",
            "} = ",saidaz2()$Statistic,"$$")))})
  output$m2conclpz2=renderText(paste("O valor p é de",saidaz2()$`p-value`, 
                                     "e o nível de significância é de",input$m2sigz21,
                                     ". Assim, conclue-se que:",
                                     ifelse(saidaz2()$`p-value`>input$m2sigz21,
                                            "Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$m2ztabe2 = renderText({paste("Com o nível de significância é ",input$m2sigz21,
                                      "o valor crítico para um teste de hipótese",
                                      ifelse(input$m2h1z1=="≠","bicaudal","unilateral"),
                                      "é Ztab = ",saidaz2()$t_value)
    paste("O limite crítico é",
          ifelse(input$m2h1z1=="≠",paste(-round(saidaz2()$t_value,4)," e ",
                                         round(saidaz2()$t_value,4)),
                 ifelse(input$m2h1z1==">",
                        paste(round(saidaz2()$t_value,4)),
                        paste(round(saidaz2()$t_value,4)))))})
  
  #========================================================================================================
  # teste z
  #========================================================================================================
  saidaz=eventReactive(input$goz2,{
    if(input$h1z2=="≠"){alternativa="greater"}
    if(input$h1z2=="<"){alternativa="esquerda"}
    if(input$h1z2==">"){alternativa="direita"}
    testez(x1 = as.numeric(input$mediasz21),
           x2 = as.numeric(input$mediasz22),
           var1=as.numeric(input$sdz21),
           var2=as.numeric(input$sdz22),
           n1=as.numeric(input$nz21),
           n2=as.numeric(input$nz22),
           alpha = as.numeric(input$sigz21),
           alternative = alternativa)})
  output$testez2=renderPlot(saidaz()$normal_plot)
  output$hip1=renderUI(withMathJax(helpText(paste("$$H_0: \\mu_1",input$h0z2,"\\mu_2$$"))))
  output$hip2=renderUI(withMathJax(helpText(paste("$$H_1: \\mu_1",input$h1z2,"\\mu_2$$"))))
  output$contavalorzz=renderUI({
    withMathJax(helpText(
      paste("$$Z_{cal}=\\frac{",input$mediasz21,"-",input$mediasz22,"}{\\sqrt{","\\frac{",input$sdz21,"^2}{",input$nz21,"}+",
            "\\frac{",input$sdz22,"^2}{",input$nz21,"}}} = ",saidaz()$Statistic,"$$")))})
  output$conclpz2=renderText(paste("O valor p é de",saidaz()$`p-value`, "e o nível de significância é de",input$sigz21,". Assim, conclue-se que:",
                                   ifelse(saidaz()$`p-value`>input$sigz21,"Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$ztabe2 = renderText({paste("Com o nível de significância é ",input$sigz21, "o valor crítico para um teste de hipótese",
                                    ifelse(input$h1z2=="≠","bicaudal","unilateral"),"é Ztab = ",saidaz()$Z_standard)
    paste("O limite crítico é",
          ifelse(input$h1z2=="≠",paste(-round(saidaz()$Z_standard,4)," e ",round(saidaz()$Z_standard,4)),
                 ifelse(input$h1z2==">",
                        paste(round(saidaz()$Z_standard,4)),
                        paste(round(saidaz()$Z_standard,4)))))})
  
  #========================================================================================================
  # teste z para prop
  #========================================================================================================
  saidapp=eventReactive(input$pgoz2,{
    if(input$h1zp2=="≠"){alternativa="greater"}
    if(input$h1zp2=="<"){alternativa="esquerda"}
    if(input$h1zp2==">"){alternativa="direita"}
    testezp(x1 = as.numeric(input$pz21),
            x2 = as.numeric(input$pz22),
            n1=as.numeric(input$pnz21),
            n2=as.numeric(input$pnz22),
            alpha = as.numeric(input$psigz21),
            alternative = alternativa)})
  output$ptestez2=renderPlot(saidapp()$normal_plot)
  output$hipp1=renderUI(withMathJax(helpText(paste("$$H_0: p_1",input$h0zp2,"p_2$$"))))
  output$hipp2=renderUI(withMathJax(helpText(paste("$$H_1: p_1",input$h1zp2,"p_2$$"))))
  output$f1p=renderUI({withMathJax(helpText(paste("$$f1=\\frac{",input$pz21,"}{",input$pnz21+input$pnz22,"} = ",saidapp()$f1,"$$")))})
  output$f2p=renderUI({withMathJax(helpText(paste("$$f2=\\frac{",input$pz22,"}{",input$pnz21+input$pnz22,"} = ",saidapp()$f2,"$$")))})
  output$pp=renderUI({withMathJax(helpText(paste("$$\\hat{p}=\\frac{",input$pz21+input$pz22,"}{",input$pnz21+input$pnz22,"} = ",saidapp()$p,"$$")))})
  
  output$pcontaz=renderUI({withMathJax(helpText(
    paste("$$Z=\\frac{",saidapp()$f1,"-",
          saidapp()$f2,"}{",saidapp()$p,
          "(1-",saidapp()$p,
          ")\\sqrt{\\frac{1}{",input$pnz21,"}+\\frac{1}{",input$pnz22,"}}} =",saidapp()$Statistic,"$$")))})
  output$pconcl2=renderText(paste("O valor p é de",saidapp()$`p-value`, "e o nível de significância é de",input$psigz21,". Assim, conclue-se que:",
                                  ifelse(saidapp()$`p-value`>input$psigz21,"Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$pztabe2 = renderText({paste("Com o nível de significância é ",input$psigz21, "o valor crítico para um teste de hipótese",
                                     ifelse(input$h1zp2=="≠","bicaudal","unilateral"),"é Ztab = ",saidapp()$Z_standard)
    paste("O limite crítico é",
          ifelse(input$h1zp2=="≠",paste(-round(saidapp()$Z_standard,4),round(saidapp()$Z_standard,4)),
                 ifelse(input$h1zp2==">",
                        paste(round(saidapp()$Z_standard,4)),
                        paste(round(saidapp()$Z_standard,4)))))})
  
  #========================================================================================================
  # teste F
  #========================================================================================================
  saidafh=eventReactive(input$govar2,{
    testefhomog(s1 = as.numeric(input$var21),
                s2 = as.numeric(input$var22),
                n1=as.numeric(input$nvar21),
                n2=as.numeric(input$nvar22),
                sig = as.numeric(input$sigvar))})
  output$contavalorvar=renderUI({
    withMathJax(helpText(
      paste("$$F_{cal}=\\frac{",saidafh()$nume,"}{",saidafh()$den,"} = ",saidafh()$fmax,"$$")))})
  output$GL1=renderUI({withMathJax(helpText(paste("$$GL_{numerador} = ",(saidafh()$gl1+1),"- 1 = ",saidafh()$gl1,"$$")))})
  output$GL2=renderUI({withMathJax(helpText(paste("$$GL_{denominador} = ",(saidafh()$gl2+1),"- 1 = ",saidafh()$gl2,"$$")))})
  output$conclpvar=renderText(paste("O valor p é de",saidafh()$p, "e o nível de significância é de",input$sigvar,". Assim, conclue-se que:",
                                    ifelse(saidafh()$p>input$sigvar,"Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$vartabe = renderText({paste("Com o nível de significância é ",input$sigvar, "o valor crítico para um teste de hipótese bilateral",
                                     "é Ztab = ",saidafh()$fat1,"a",saidafh()$fat2)
    paste("O limite crítico inferior é ",
          round(saidafh()$ft2,4), "e o limite crítico superior é",round(saidafh()$ft1,4))})
  
  #========================================================================================================
  # teste t
  #========================================================================================================
  saidat=eventReactive(input$got2,{
    if(input$h1t2=="≠"){alternativa="greater"}
    if(input$h1t2=="<"){alternativa="esquerda"}
    if(input$h1t2==">"){alternativa="direita"}
    testet_vari(x1 = as.numeric(input$mediast21),
                x2 = as.numeric(input$mediast22),
                s1=as.numeric(input$sdt21),
                s2=as.numeric(input$sdt22),
                n1=as.numeric(input$nt21),
                n2=as.numeric(input$nt22),
                alpha = as.numeric(input$sigt21),
                alternative = alternativa)})
  output$testetg2=renderPlot(saidat()$t_plot)
  output$hipt1=renderUI(withMathJax(helpText(paste("$$H_0: \\mu_1",input$h0t2,"\\mu_2$$"))))
  output$hipt2=renderUI(withMathJax(helpText(paste("$$H_1: \\mu_1",input$h1t2,"\\mu_2$$"))))
  output$scvalortt=renderUI({
    withMathJax(helpText(
      paste("$$S_c^2=\\frac{(",input$nt21,"-1)",input$sdt21,"^2+(",input$nt22,"-1)",input$sdt22,"^2}{",input$nt21+input$nt21-2,"} =",saidat()$var_c,"$$")))})
  output$contavalortt=renderUI({
    withMathJax(helpText(
      paste("$$t=\\frac{",input$mediast21," - ",input$mediast22,"}{",
            sqrt(saidat()$var_c),"\\sqrt{\\frac{1}{",input$nt21,"}+\\frac{1}{",input$nt22,"}}} = ",saidat()$Statistic,"$$")))})
  
  output$conclpt2=renderText(paste("O valor p é de",saidat()$`p-value`, "e o nível de significância é de",input$sigt21,". Assim, conclue-se que:",
                                   ifelse(saidat()$`p-value`>input$sigt21,"Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$ttabe2 = renderText({paste("Com o nível de significância é ",input$sigt21, "o valor crítico para um teste de hipótese",
                                    ifelse(input$h1t2=="≠","bicaudal","unilateral"),"é Ztab = ",saidat()$t_value)
    paste("O limite crítico é ",
          ifelse(input$h1t2=="≠",paste(-round(saidat()$t_value,4)," e ",round(saidat()$t_value,4)),
                 ifelse(input$h1t2==">",
                        paste(round(saidat()$t_value,4)),
                        paste(round(saidat()$t_value,4)))))})
  
  #========================================================================================================
  # teste t 2
  #========================================================================================================
  saidatt=eventReactive(input$gott2,{
    if(input$h1tt2=="≠"){alternativa="greater"}
    if(input$h1tt2=="<"){alternativa="esquerda"}
    if(input$h1tt2==">"){alternativa="direita"}
    testet_vard(x1 = as.numeric(input$mediastt21),
                x2 = as.numeric(input$mediastt22),
                s1=as.numeric(input$sdtt21),
                s2=as.numeric(input$sdtt22),
                n1=as.numeric(input$ntt21),
                n2=as.numeric(input$ntt22),
                alpha = as.numeric(input$sigtt21),
                alternative = alternativa)})
  output$testettg2=renderPlot(saidatt()$t_plot)
  output$hiptt1=renderUI(withMathJax(helpText(paste("$$H_0: \\mu_1",input$h0tt2,"\\mu_2$$"))))
  output$hiptt2=renderUI(withMathJax(helpText(paste("$$H_1: \\mu_1",input$h1tt2,"\\mu_2$$"))))
  output$contavalorttt=renderUI({
    withMathJax(helpText(
      paste("$$t_{cal}=\\frac{",input$mediastt21,"-",input$mediastt22,"}{\\sqrt{","\\frac{",input$sdtt21,"^2}{",input$ntt21,"}+",
            "\\frac{",input$sdtt22,"^2}{",input$ntt21,"}}} = ",saidatt()$Statistic,"$$")))})
  output$vvalorttt=renderUI({
    withMathJax(helpText(
      paste("$$v = \\frac{(\\frac{",input$sdtt21,
            "^2}{",input$ntt21,
            "}+\\frac{",input$sdtt22,
            "^2}{",input$ntt22,
            "})^2}{\\frac{(",input$sdtt21,
            "^2/",input$ntt21,
            ")^2}{",input$ntt21,
            "-1}+\\frac{(",input$sdtt22,
            "^2/",input$ntt22,
            ")^2}{",input$ntt22,"-1}} = ",saidatt()$v,"$$")))})
  output$conclptt2=renderText(paste("O valor p é de",saidatt()$`p-value`, "e o nível de significância é de",input$sigtt21,". Assim, conclue-se que:",
                                    ifelse(saidatt()$`p-value`>input$sigtt21,"Não rejeita a hipotese nula","Rejeita-se a hipótese nula")))
  output$tttabe2 = renderText({paste("Com o nível de significância é ",input$sigtt21, "o valor crítico para um teste de hipótese",
                                     ifelse(input$h1tt2=="≠","bicaudal","unilateral"),"é Ztab = ",saidatt()$t_value)
    paste("O limite crítico é",
          ifelse(input$h1tt2=="≠",paste(-round(saidatt()$t_value,4)," e ",round(saidatt()$t_value,4)),
                 ifelse(input$h1tt2==">",
                        paste(round(saidatt()$t_value,4)),
                        paste(round(saidatt()$t_value,4)))))})
  output$tabelatestet=renderDataTable({round(tabelat,3)}, options = list(pageLength = 50))
  output$tabelatestez=renderDataTable({round(tabelaz,4)}, options = list(pageLength = 40))
  output$tabelatestef=renderDataTable({tabelaf(sig = input$sigF)}, options = list(pageLength = 37))
  
  #========================================================================================================
  # Anova DIC
  #========================================================================================================
  z <- reactive({
    infile <- input$uploadedcsvaov
    if (is.null(infile)){return(cat("Adicione seu conjunto de dados"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhaaov))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvaov)}})
  observe({
    vchoices <- names(z())
    updateSelectInput(session, "tratA", choices = vchoices,selected = names(z())[1])
    updateSelectInput(session, "respA", choices = vchoices,selected = names(z())[2])})
  dadosaov=eventReactive(input$goaov,{
    infile <- input$uploadedcsvaov
    if (is.null(infile)){return(cat("Selecione a coluna de tratamento e resposta"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhaaov))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvaov)}
    trat=tbl[,input$tratA]
    resp=tbl[,input$respA]
    nlevels=length(levels(as.factor(trat)))
    trat=factor(trat,unique(trat))
    soma=tapply(resp,trat,sum)
    n=length(resp)/length(levels(trat))
    CC=round(sum(resp)^2/length(resp),5)
    sqtrat=sum(soma^2)/n-CC
    sqtotal=sum(resp^2)-CC
    sqres=sqtotal-sqtrat
    gltotal=length(resp)
    gltrat=nlevels
    list("trat"=trat,
         "resp"=resp,
         "nlevels"=nlevels,
         soma=soma,
         n=n,
         CC=CC,sqtrat=sqtrat,
         sqtotal=sqtotal,
         sqres=sqres,
         gltotal=gltotal)})
  output$correcao1=renderUI({
    resp=dadosaov()$resp
    n=length(resp)
    p1=paste(resp[1:3],rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),resp[(n-2):n],collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$ C = \\frac{(",
                               final,
                               ")^2}{",
                               length(dadosaov()$resp),"} = ",
                               round(sum(dadosaov()$resp)^2/length(dadosaov()$resp),5),"$$")))})
  output$gltotal=renderUI({
    withMathJax(helpText("$$ GL_{total} = ",dadosaov()$gltotal,"- 1 = ",dadosaov()$gltotal-1,"$$"))})
  output$gltrat=renderUI({
    withMathJax(helpText("$$ GL_{trat} = ",dadosaov()$nlevels,"- 1 = ",dadosaov()$nlevels-1,"$$"))})
  output$glresiduo=renderUI({
    withMathJax(helpText("$$ GL_{res} = ",dadosaov()$gltotal-1,"-",dadosaov()$nlevels-1, " = ",dadosaov()$gltotal-dadosaov()$nlevels,#-2,
                         "$$"))})
  output$sqt=renderUI({
    n=dadosaov()$gltotal
    p1=paste(paste(dadosaov()$resp[1:3],"^2"),rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),paste(dadosaov()$resp[(n-2):n],"^2"),collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$ SQT = ",
                               final,
                               "-",dadosaov()$CC,
                               " = ",
                               round(dadosaov()$sqtotal,5),"$$")))})
  output$sqtrat=renderUI({
    withMathJax(helpText(paste("$$ SQ_{trat} = \\frac{",
                               paste(paste(dadosaov()$soma,"^2"),c(rep(" + ",length(dadosaov()$soma)-1),""),
                                     collapse = "",sep=""),"}{",dadosaov()$n,"}",
                               "-",dadosaov()$CC,
                               " = ",
                               round(dadosaov()$sqtrat,5),"$$")))})
  output$sqres=renderUI({
    withMathJax(helpText(paste("$$ SQ_{res} = ",round(dadosaov()$sqtotal,4), " - ", 
                               round(dadosaov()$sqtrat,4), " = ",
                               round(dadosaov()$sqtotal-dadosaov()$sqtrat,4),"$$")))})
  output$qmtrat=renderUI({
    withMathJax(helpText(paste("$$ QM_{trat} = \\frac{",round(dadosaov()$sqtrat,4), "}{",(dadosaov()$nlevels-1) ,"} = ",
                               round(dadosaov()$sqtrat/(dadosaov()$nlevels-1),4),"$$")))})
  output$qmres=renderUI({
    gltotal=length(dadosaov()$resp)-1
    gltrat=dadosaov()$nlevels-1
    withMathJax(helpText(paste("$$ QM_{res} = \\frac{",round(dadosaov()$sqres,4), "}{",gltotal-gltrat ,"} = ",
                               round(dadosaov()$sqres/(gltotal-gltrat),4),"$$")))})
  output$fcal=renderUI({
    gltotal=length(dadosaov()$resp)-1
    gltrat=dadosaov()$nlevels-1
    glres=gltotal-gltrat
    withMathJax(helpText(paste("$$ F_{cal} = \\frac{",round(dadosaov()$sqtrat/gltrat,4), "}{",round(dadosaov()$sqres/glres,4),"} = ",
                               round((dadosaov()$sqtrat/gltrat)/(dadosaov()$sqres/glres),4),"$$")))})
  output$quadroanovadic = shiny::renderDataTable({
    ano=aov(dadosaov()$resp~as.factor(dadosaov()$trat))
    ano1=data.frame(anova(ano))
    ano1=cbind(trat=c("Tratamento","Resíduo"),ano1)
    colnames(ano1)=c("FV","GL","SQ","QM","Fc","p-valor")
    ano1=rbind(ano1,c("Total",sum(ano1$GL),sum(ano1$SQ),sum(ano1$QM),NA,NA))
    ano1})
  output$dmstukeydic=renderUI({
    r=dadosaov()$n
    gltotal=length(dadosaov()$resp)-1
    gltrat=dadosaov()$nlevels-1
    glres=gltotal-gltrat
    s=dadosaov()$sqres/glres
    qtuk = qtukey(1-input$alphatukeydic,nmeans = dadosaov()$nlevels,glres) 
    dms=qtuk*sqrt(s/r)
    withMathJax(helpText(paste("$$ \\Delta = ",qtuk,"\\times \\sqrt{\\frac{",s,"}{",r,"}} = ",dms,"$$")))})
  output$tukeytabela=shiny::renderDataTable({
    r=dadosaov()$n
    gltotal=length(dadosaov()$resp)-1
    gltrat=dadosaov()$nlevels-1
    glres=gltotal-gltrat
    s=dadosaov()$sqres/glres
    a=agricolae::HSD.test(dadosaov()$resp,dadosaov()$trat,DFerror = glres,MSerror = s,alpha = input$alphatukeydic)
    nomes=rownames(a$groups)
    final=cbind("Trat"=nomes,a$groups)
    colnames(final)=c("Trat","Médias","Tukey")
    final})
  
  #========================================================================================================
  # Anova DBC
  #========================================================================================================
  z1 <- reactive({
    infile <- input$uploadedcsvaov1
    if (is.null(infile)){return(cat("Adicione seu conjunto de dados"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhaaov1))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvaov1)}})
  observe({
    vchoices <- names(z1())
    updateSelectInput(session, "tratB", choices = vchoices,selected = names(z1())[1])
    updateSelectInput(session, "blocoB", choices = vchoices,selected = names(z1())[2])
    updateSelectInput(session, "respB", choices = vchoices,selected = names(z1())[3])})
  dadosaov1=eventReactive(input$goaov1,{
    infile <- input$uploadedcsvaov1
    if (is.null(infile)){return(cat("Selecione a coluna de tratamento, bloco e resposta"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhaaov1))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvaov1)}
    trat=tbl[,input$tratB]
    bloco=tbl[,input$blocoB]
    resp=tbl[,input$respB]
    nlevels=length(levels(as.factor(trat)))
    nlevelsbloco=length(levels(as.factor(bloco)))
    trat=factor(trat,unique(trat))
    bloco=factor(bloco,unique(bloco))
    somatrat=tapply(resp,trat,sum)
    somabloco=tapply(resp,bloco,sum)
    n=length(resp)/length(levels(trat))
    CC=round(sum(resp)^2/length(resp),5)
    sqtrat=sum(somatrat^2)/nlevelsbloco-CC
    sqbloco=sum(somabloco^2)/nlevels-CC
    sqtotal=sum(resp^2)-CC
    sqres=sqtotal-sqtrat-sqbloco
    gltotal=length(resp)
    gltrat=nlevels
    glbloco=nlevelsbloco
    list("trat"=trat,
         "bloco"=bloco,
         "resp"=resp,
         "nlevels"=nlevels,
         "nlevelsbloco"=nlevelsbloco,
         "somatrat"=somatrat,
         "somabloco"=somabloco,
         n=n,
         CC=CC,
         sqtrat=sqtrat,
         sqbloco=sqbloco,
         sqtotal=sqtotal,
         sqres=sqres,
         gltotal=gltotal)})
  output$correcao2=renderUI({
    resp=dadosaov1()$resp
    n=length(resp)
    p1=paste(resp[1:3],rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),resp[(n-2):n],collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$ C = \\frac{(",
                               final,
                               ")^2}{",
                               length(dadosaov1()$resp),"} = ",
                               round(sum(dadosaov1()$resp)^2/length(dadosaov1()$resp),5),"$$")))})
  output$gltotal1=renderUI({
    withMathJax(helpText("$$ GL_{total} = ",dadosaov1()$gltotal,"- 1 = ",dadosaov1()$gltotal-1,"$$"))})
  output$gltrat1=renderUI({
    withMathJax(helpText("$$ GL_{trat} = ",dadosaov1()$nlevels,"- 1 = ",dadosaov1()$nlevels-1,"$$"))})
  output$glbloc1=renderUI({
    withMathJax(helpText("$$ GL_{bloco} = ",dadosaov1()$nlevelsbloco,"- 1 = ",dadosaov1()$nlevelsbloco-1,"$$"))})
  output$glresiduo1=renderUI({
    withMathJax(helpText("$$ GL_{res} = ",
                         dadosaov1()$gltotal-1,"-",
                         dadosaov1()$nlevels-1, "-",
                         dadosaov1()$nlevelsbloco-1," = ",
                         (dadosaov1()$gltotal-1)-(dadosaov1()$nlevels-1)-
                           (dadosaov1()$nlevelsbloco-1),#-2,
                         "$$"))})
  output$sqt1=renderUI({
    n=dadosaov1()$gltotal
    p1=paste(paste(dadosaov1()$resp[1:3],"^2"),rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),paste(dadosaov1()$resp[(n-2):n],"^2"),collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$ SQT = ",
                               final,
                               "-",dadosaov1()$CC,
                               " = ",
                               round(dadosaov1()$sqtotal,5),"$$")))})
  output$sqbloc1=renderUI({
    withMathJax(helpText(paste("$$ SQ_{bloco} = \\frac{",
                               paste(paste(dadosaov1()$somabloco,"^2"),
                                     c(rep(" + ",length(dadosaov1()$somabloco)-1),""),
                                     collapse = "",sep=""),"}{",dadosaov1()$nlevels,"}",
                               "-",dadosaov1()$CC,
                               " = ",
                               round(dadosaov1()$sqbloco,5),"$$")))})
  output$sqtrat1=renderUI({
    withMathJax(helpText(paste("$$ SQ_{trat} = \\frac{",
                               paste(paste(dadosaov1()$somatrat,"^2"),
                                     c(rep(" + ",length(dadosaov1()$somatrat)-1),""),
                                     collapse = "",sep=""),"}{",dadosaov1()$nlevelsbloco,"}",
                               "-",dadosaov1()$CC,
                               " = ",
                               round(dadosaov1()$sqtrat,5),"$$")))})
  output$sqres1=renderUI({
    withMathJax(helpText(paste("$$ SQ_{res} = ",round(dadosaov1()$sqtotal,4), " - ",
                               round(dadosaov1()$sqtrat,4)," - ",
                               round(dadosaov1()$sqbloco,4), " = ",
                               round(dadosaov1()$sqtotal-dadosaov1()$sqtrat-dadosaov1()$sqbloco,4),"$$")))})
  output$qmtrat1=renderUI({
    withMathJax(helpText(paste("$$ QM_{trat} = \\frac{",
                               round(dadosaov1()$sqtrat,4), "}{",
                               (dadosaov1()$nlevels-1) ,"} = ",
                               round(dadosaov1()$sqtrat/(dadosaov1()$nlevels-1),4),"$$")))})
  output$qmbloc1=renderUI({
    withMathJax(helpText(paste("$$ QM_{bloco} = \\frac{",
                               round(dadosaov1()$sqbloco,4), "}{",(dadosaov1()$nlevelsbloco-1) ,"} = ",
                               round(dadosaov1()$sqbloco/(dadosaov1()$nlevelsbloco-1),4),"$$")))})
  output$qmres1=renderUI({
    gltotal=length(dadosaov1()$resp)-1
    gltrat=dadosaov1()$nlevels-1
    glbloco=dadosaov1()$nlevelsbloco-1
    withMathJax(helpText(paste("$$ QM_{res} = \\frac{",round(dadosaov1()$sqres,4), "}{",
                               gltotal-gltrat-glbloco ,"} = ",
                               round(dadosaov1()$sqres/(gltotal-gltrat-glbloco),4),"$$")))})
  output$fcal1=renderUI({
    gltotal=length(dadosaov1()$resp)-1
    gltrat=dadosaov1()$nlevels-1
    glbloco=dadosaov1()$nlevelsbloco-1
    glres=gltotal-gltrat-glbloco
    withMathJax(helpText(paste("$$ F_{trat} = \\frac{",
                               round(dadosaov1()$sqtrat/gltrat,4), "}{",
                               round(dadosaov1()$sqres/glres,4),"} = ",
                               round((dadosaov1()$sqtrat/gltrat)/
                                       (dadosaov1()$sqres/glres),4),"$$")))})
  output$fcalbloco1=renderUI({
    gltotal=length(dadosaov1()$resp)-1
    gltrat=dadosaov1()$nlevels-1
    glbloco=dadosaov1()$nlevelsbloco-1
    glres=gltotal-gltrat-glbloco
    withMathJax(helpText(paste("$$ F_{bloco} = \\frac{",
                               round(dadosaov1()$sqbloco/glbloco,4), "}{",
                               round(dadosaov1()$sqres/glres,4),"} = ",
                               round((dadosaov1()$sqbloco/glbloco)/
                                       (dadosaov1()$sqres/glres),4),"$$")))})
  output$quadroanovadbc = shiny::renderDataTable({
    ano=aov(dadosaov1()$resp~as.factor(dadosaov1()$trat)+as.factor(dadosaov1()$bloco))
    ano1=data.frame(anova(ano))
    ano1=cbind(trat=c("Tratamento","Bloco","Resíduo"),ano1)
    colnames(ano1)=c("FV","GL","SQ","QM","Fc","p-valor")
    ano1=rbind(ano1,c("Total",sum(ano1$GL),sum(ano1$SQ),sum(ano1$QM),NA,NA))
    ano1})
  output$dmstukeydbc=renderUI({
    r=dadosaov1()$nlevelsbloco
    gltotal=length(dadosaov1()$resp)-1
    gltrat=dadosaov1()$nlevels-1
    glbloco=dadosaov1()$nlevelsbloco-1
    glres=gltotal-gltrat-glbloco
    s=dadosaov1()$sqres/glres
    qtuk = qtukey(1-input$alphatukeydbc,nmeans = dadosaov1()$nlevels,glres) 
    dms=qtuk*sqrt(s/r)
    withMathJax(helpText(paste("$$ \\Delta = ",qtuk,"\\times \\sqrt{\\frac{",s,"}{",r,"}} = ",dms,"$$")))})
  output$tukeytabeladbc=shiny::renderDataTable({
    r=dadosaov1()$nlevelsbloco
    gltotal=length(dadosaov1()$resp)-1
    gltrat=dadosaov1()$nlevels-1
    glbloco=dadosaov1()$nlevelsbloco-1
    glres=gltotal-gltrat-glbloco
    s=dadosaov1()$sqres/glres
    a=agricolae::HSD.test(dadosaov1()$resp,dadosaov1()$trat,DFerror = glres,MSerror = s,alpha = input$alphatukeydbc)
    nomes=rownames(a$groups)
    final=cbind("Trat"=nomes,a$groups)
    colnames(final)=c("Trat","Médias","Tukey")
    final})
  
  #========================================================================================================
  # Correlacao
  #========================================================================================================
  corre1 <- reactive({
    infile <- input$uploadedcsvcor
    if (is.null(infile)){return(cat("Adicione seu conjunto de dados"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhacor))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvcor)}})
  observe({
    vchoices <- names(corre1())
    updateSelectInput(session, "X", choices = vchoices,selected = names(z1())[1])
    updateSelectInput(session, "Y", choices = vchoices,selected = names(z1())[2])})
  dadoscorre=eventReactive(input$gocor,{
    infile <- input$uploadedcsvcor
    if (is.null(infile)){return(cat("Selecione a coluna das variáveis"))}
    ext=tools::file_ext(infile)
    if(ext[1]=="xlsx"| ext[1]=="xls"){
      tbl=data.frame(readxl::read_excel(infile$datapath,sheet = input$planilhacor))}
    else{tbl=read.csv(infile$datapath, header = TRUE, sep = input$deccsvcor)}
    trat=tbl[,input$X]
    resp=tbl[,input$Y]
    list("trat"=trat,
         "resp"=resp)})
  output$mediaX=renderUI({
    n=length(dadoscorre()$trat)
    p1=paste(dadoscorre()$trat[1:3],rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),dadoscorre()$trat[(n-2):n],collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$\\bar{x} = \\frac{",final,"}{",length(dadoscorre()$trat),"} = ",
                               mean(dadoscorre()$trat),"$$")))})
  output$mediaY=renderUI({
    n=length(dadoscorre()$resp)
    p1=paste(dadoscorre()$resp[1:3],rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),dadoscorre()$resp[(n-2):n],collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$\\bar{y} = \\frac{",final,"}{",length(dadoscorre()$resp),"} = ",
                               mean(dadoscorre()$resp),"$$")))})
  output$desvioX=renderUI({
    n=length(dadoscorre()$trat)
    p1=paste(paste(dadoscorre()$trat[1:3],"^2"),rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),paste(dadoscorre()$trat[(n-2):n],"^2"),collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$S_{x} = \\sqrt{\\frac{",final,"}{",
                               length(dadoscorre()$trat),"- 1}} = ",
                               sd(dadoscorre()$trat),"$$")))})
  output$desvioY=renderUI({
    n=length(dadoscorre()$resp)
    p1=paste(paste(dadoscorre()$resp[1:3],"^2"),rep(" + ",3),collapse = "",sep = "")
    p2=paste(rep(" + ",3),paste(dadoscorre()$resp[(n-2):n],"^2"),collapse = "",sep = "")
    final=paste(p1,"\\cdots",p2,collapse = "")
    withMathJax(helpText(paste("$$S_{y} = \\sqrt{\\frac{",final,"}{",
                               length(dadoscorre()$resp)," - 1}} = ",
                               sd(dadoscorre()$resp),"$$")))})
  output$pearson=renderUI({
    mx=mean(dadoscorre()$trat)
    my=mean(dadoscorre()$resp)
    sdx=sd(dadoscorre()$trat)
    sdy=sd(dadoscorre()$resp)
    resp=dadoscorre()$resp
    trat=dadoscorre()$trat
    p1=paste(paste("(",paste(trat[1:3],rep("-",3),mx),")"," (",paste(resp[1:3],rep("-",3),my),") +"),collapse = "")
    n=length(resp)
    nfim=(n-2)
    p2=paste(paste("(",paste(resp[nfim:n],rep("-",3),mx),")"," (",paste(resp[nfim:n],rep("-",3),my),") +"),collapse = "")
    final=paste(p1,"\\cdots",p2,collapse = "")  
    withMathJax(helpText(paste("$$r = \\frac{",final,"}{(",length(dadoscorre()$trat),"-1)",sdx,"\\times",sdy,"} = ",
                               cor(trat,resp),"$$")))})
}

shinyApp(ui, server)
