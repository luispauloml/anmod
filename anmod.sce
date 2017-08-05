//====================================
// 
// DISTRIBUTED UNDER GPLv3
// 
//====================================
// 
// Luis Paulo LIMA         mar/2016
// Engenharia Mecânica
// Universidade Federal de Mato Grosso, Brasil
// <luispauloml[arroba]gmail.com>
// 
// Disponível em: <http://pastebin.com/TUmmg18U>
//
//====================================
// Análise modal de sist/emas vibratórios com amortecimento histerético
// proprocional com estudo teório e experimental, usando o método da 
// amplitude dos picos para para extração de parâmetros modais.
// 
// [Entrada:]
// - Modelo espacial:
//         < Rota teórica >
//     -- N            número de graus de liberdade
//     -- M            matriz de massa
//     -- K            matriz de rigidez
//     -- betah        proporção de amortecimento relativo à rigidez
//         < Rota experimental >
//     -- n            índice da linha para qual se reduz a 
//                     resposta do sistema
//     
// [Saída:]
//         < Rota teórica >
// - Amortecimento e modelo espacial generalizado:
//     -- D            matriz de amortecimento
//     -- m_r_diag     matriz de massa generalizada
//     -- k_r_diag     matriz de rigidez generalizada
//     -- d_r_diag     matriz de amortecimento generalizado
// - Modelo modal:
//     -- w_n          vetor de frequências natruais
//     -- etah_n       vetor de fatores de amortecimento modal
//     -- PSI          forma dos modos normais
//     -- PHI          forma dos modos normais normalizados pela massa
// - Modelo de resposta
//     -- w           vetor de frequências para cálculo da resposta
//     -- alfa1       matriz de expressões simbólicas das funções 
//                    receptânicas
//     -- H           hipermatriz com valores das receptâncias para
//                    cada frequência
//                     < H(j,:,k) é o vetor que contem os valores 
//                     calculados de alfa1(j,k) >
// 
//         < Rota experimental >
// - Modelo modal experimental
//     -- w_n_med      vetor das frequências naturais medidas 
//                     graficamente
//     -- etah_n_med   vetor dos fatores de amortecimento modal 
//                     medidos graficamente
//     -- psi_graf     forma dos modos normais experimentais
//     -- phi_graf     forma dos modos normais normalizados pela 
//                     massa experimentais
// - Modelo espacial experimental
//     -- M_exp        matriz de massa experimental
//     -- K_exp        matriz de rigidez experimetal
//     -- betah_exp    proporção de amortecimento experimental
//     -- D_exp        matriz de amortecimento experimental
//     -- m_r_exp      matriz de massa generalizada experimental
//     -- k_r_exp      matriz de rigidez generalizada experimental
//     -- d_r_exp      matriz de amortecimento generalizado 
//                     experimental
// - Modelo de resposta
//     -- alfa2        matriz de expressões simbólicas das funções 
//                     receptânicas
//     -- H2           hipermatriz com valores das receptâncias para 
//                     cada frequência
//                     < H2(j,:,k) é o vetor que contem os valores 
//                     calculados de alfa2(j,k) >
// 
//         < Gráficos >
// -- Janela 10    diagramas de Bode para alfa1(1,1) a alfa1(1,3)
//                 picos destacados
//                 pontos de meia potência destacados
// -- Janela 20    diagramas de Bode para alfa1(1,3) a alfa1(1,6)
//                 picos destacados
//                 pontos de meia potência destacados
// -- Janela 11    diagrama de Nyquist para alfa1(1,1) a alfa1(1,3)
// -- Janela 30    comparação da curva original alfa11(1,4) e sua
//                 versão regenerada alfa2(1,4)
//======FIM=====



//============ [ R O T A    T E Ó R I C A ] ============

//1. DADOS INICIAIS
N=6;    //número de graus de liberdade [admens.]
massas=[7 7 4 3 6 8];   //massas [kg]

k1=1.0d5; k2=1.0d5;  k3=4.0d5;  k4=5.0d5;
k5=7.0d5; k6=2.0d5;  k7=8.0d5;  k8=3.0d5;
k9=6.0d5; k10=3.0d5; k11=5.0d5;

molas=[k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11] //matriz das molas [N/m]

betah=1.0d-3;            //proporção de amortecimento relativo
                         //à rigidez [admens.]
    
//2. MODELO ESPACIAL ANALÍTICO
M=[]; //Matriz de massas [kg]
K=[]; //matriz de rigidez [N/m]
D=[]; //matriz de amortecimento [N/m]


    M=diag(massas(1:N));        //matriz de massas
    
    K=...                       //matriz de rigidez
    [k1+k3+k6  0  -k3  -k6  0  0;...
        0  k2+k4+k5  -k4  0  -k5  0;...
       -k3  -k4  k3+k4+k7+k8+k9  -k8  -k7  -k9;...
       -k6  0  -k8  k6+k8+k11  0  -k11;...
        0  -k5  -k7  0  k5+k7+k10  -k10;...
        0  0  -k9  -k11  -k10  k9+k10+k11];
    
    D=betah*K;                  //matriz de amortecimento

//3. MODELO MODAL ANALÍTICO
//3.1 Autoproblema
PSI=[];                 //Modos de vibração [m]
w_n=[];                 //Frequências naturais [rad/s]
etah_n=[];              //Fatores de amortecimento [admens.]
lambda_n2=[];           //Matriz de autovalores ao quadrado 
                        //[lambda2=w_n2(1+i*etah)]

    [al,be,Q]=spec(K+%i*D,M);   //Resolver autoproblema
    autov=al./be;               //Autovalores desordenados
    
    clear al be
    
//3.1.1 Ordenar autovetores e autovalores
    ReAutov=real(autov);
    [ReAutov,ind]=gsort(ReAutov,"g","i");
    lambda_n2=autov(ind);   //Autovalores

    w_n=sqrt(ReAutov); w_n=real(w_n);   //Frequências naturais
    etah_n=imag(lambda_n2)./ReAutov;    //Fatores de amortecimento
    
    Q1=zeros(Q);        //Autovetores sem partes imaginárias
    tol=1.0d-3;         //Tolerância da parte imaginária
    for r=1:size(Q,"r") //Eliminando partes imaginárias da matriz 
                        //de autovetores
        for c=1:size(Q,"c")
            if (imag(Q(r,c))/real(Q(r,c))<tol)|...
                (abs(imag(Q(r,c))/real(Q(r,c))-betah)<tol) then
                Q1(r,c)=real(Q(r,c));
            else
                Q1(r,c)=Q(r,c);
                print(%io(2),[r c])
            end
        end
    end
    
    PSI=Q1(:,ind);          //Modos de vibrar ordenados em 
                            //relação às frequências naturais
                            
    for c=1:size(PSI,"c")   //Normalização em relação a coordenada 1
        PSI(:,c)=PSI(:,c)./real(PSI(1,c));
    end
    
    clear ind r c tol
    
//3.2 Normalizar pela massa
m_r_diag=[];    //matriz massa diagonalziada [kg]
k_r_diag=[];    //matriz rigidez diagonalizada [N/m]
d_r_diag=[];    //matriz de amortecimento diagonalizada [Nm/s]
PHI=[];         //autovetores normalizados pela massa

    m_r_diag=PSI'*M*PSI;    //massa diagonalizada
    k_r_diag=PSI'*K*PSI;    //rigidez diagonalizada
    d_r_diag=PSI'*D*PSI;    //amortecimento diagonalizado
    
    for i=1:size(m_r_diag,"r");
        m_r(i)=m_r_diag(i,i);       //eliminando elementos não-zero 
                                    //(muito pequenos)...
        k_r(i)=k_r_diag(i,i);       //garantindo autovetores normalizados
                                    //reais
        d_r(i)=d_r_diag(i,i);
    end
    
    m_r_diag=diag(m_r);
    k_r_diag=diag(k_r);
    d_r_diag=diag(d_r);
    
    PHI=PSI*sqrt(inv(m_r_diag));   //autovetores normalizados pela massa
    
    clear i
    
//4. MODELO DE RESPOSTA ANALÍTICO
//4.1 Funções de resposta em frequência (FRF)
alfa1=[];      //matriz de FRF (receptãncia)
    
    alfa1=zeros(N,N);
    
    for j=1:N
        for k=1:N
            for r=1:N
                alfa1(j,k) = alfa1(j,k) + PHI(j,r)*PHI(k,r) ./ ...
                ( w_n(r)^2 - %z^2 + %i*etah_n(r)*w_n(r)^2 );
            end
        end
    end
        
//6. CALCULANDO AS FRFs
    w=[];       //Matriz de frequências
    H=[];       //Matriz de receptância
    
    min_w=min(w_n)-50;
    if min_w<0 then
        min_w=0.001;
    end
    
    max_w=max(w_n)+100;
    
    w=linspace(min_w,-2+w_n(1),100);
    for i=1:N-1
        w=[w(1:length(w)-1) linspace(-2+w_n(i),2+w_n(i),200)];
        w=[w(1:length(w)-1) linspace(w(length(w)),-2+w_n(i+1),100)];
    end
    w=[w(1:length(w)-1) linspace(-2+w_n(N),2+w_n(N),200)];
    w=[w(1:length(w)-1) linspace(2+w_n(N),max_w,100)];
    
    H=zeros(N,length(w),N);
    
    
    for j=1:N
        for k=1:N
            H(j,:,k)=horner(alfa1(j,k),w);
        end
    end
    
    clear A j k
    
//X. EXIBIÇÃO DOS MODELOS
    disp('[ROTA TEÓRICA]');
    disp('Modelo espacial');
    print(%io(2),D,K,M);
    
    disp(' ')
    disp('Modelo modal');
    print(%io(2),PHI,PSI,etah_n',w_n');
    
    disp(' ')
    disp('Modelo de reposta (alfa_11)');
    disp(alfa1(1,1))

//X.1 Gráficos
    scf(10); clf();
    xname("FRF %d");
    subplot(211);
    plot2d('nl', w', abs([H(1,:,1)' H(1,:,2)' H(1,:,3)']),...
           [2,3,5],leg="H11@H12@H13");
    subplot(212);
    plot2d(w',[phasemag([H(1,:,1); H(1,:,2); H(1,:,3)],'m')]',[2,3,5]);
    
    scf(20); clf();
    xname("FRF %d");
    subplot(211);
    plot2d('nl', w', abs([H(1,:,4)' H(1,:,5)' H(1,:,6)']),...
            [2,3,5],leg="H14@H15@H16");
    subplot(212);
    plot2d(w',[phasemag([H(1,:,4);H(1,:,5);H(1,:,6)],'m')]',[2,3,5]);
    
    scf(11); clf();
    xname("FRF Nyquist %d");
    plot2d4(real([H(1,:,1)' H(1,:,2)' H(1,:,2)']),...
            imag(-[H(1,:,1)' H(1,:,2)' H(1,:,2)']));
    eixo=gca();
    eixo.isoview="on";



//============ [ R O T A    E X P E R I M E N T A L ] ============

function peaks=peak_detect(signal,threshold)
// Disponível em: <https://fileexchange.scilab.org/toolboxes/209000>
//
// This function detect the peaks of a signal : 
// --------------------------------------------
// For an input row vector "signal" , the function return 
// the position of the peaks of the signal.
//
// The ouput "peaks" is a row vector (size = number of peaks),
// "peaks" =[] if no peak is found.
//
// Optional argument "threshold" eliminates the peaks under
// the threshold value (noise floor). 
//
// Clipped peaks (more than 2 samples of the signal at the same value)
// are not detected.
// -------------------------------------------------------------------
//     Jean-Luc GOUDIER      11-2011
// -------------------------------------------------------------------
    
    
    [nargout,nargin] = argn(0);
    if nargin==2 then ts=threshold;
    end;
    if nargin==1 then ts=min(signal);
    end;
    
    [r c]=size(signal);
    Ct=getlanguage();
    if Ct=="fr_FR" then
         Msg="Erreur : le signal n''est pas un vecteur colonne";
    else
         Msg="Error : signal is not a row vector";
    end
    if r>1 then
        error(Msg);
    end;
    
    Lg=c-1; 
    d_s=diff(signal); 
    dd_s=[d_s(1),d_s(1,:)];               // diff first shift
    d_s=[d_s(1,:),d_s(Lg)];               // diff size correction
    ddd_s=[dd_s(1),dd_s(1,1:Lg)];         // diff second shift
    Z=d_s.*dd_s;                          // diff zeros
    
    peaks=find(((Z<0 & d_s<0)|(Z==0 & d_s<0 & ddd_s>0)) & signal>ts);

endfunction

//7. EXTRAÇÃO DE PARÂMETROS MODAIS
n=1;                //Ponto de respota para modelo reduzido
w_n_graf=[];        //Frequências naturais extraídas
w_n_med=[];         //Valor médio das frequências naturais
etah_graf=[];       //Fator de amortecimento para cada modo
etah_med=[];        //Valor médio do amortecimento modal
lim=[];             //Limite inferior para busca de picos
lambda_n2_med=[];   //Autovalores médios

lim=5.0d-7;     //segundo inspeção gráfica

//7.1 Frequências naturais
    w_n_ind=zeros(N,N);     //índice dos picos para FRF H(n,n)
    Hnn_nat=[];             //amplitude dos picos da FRF H(n,n)
    
    for j=1:N
        w_n_ind(j,:)=peak_detect(abs(H(n,:,j)),lim);    //amplitude
        w_n_graf(j,:)=w(w_n_ind(j,:));                  //freq natural
        
        Hnn_nat(j,:)=H(n,w_n_ind(j,:),j);
    end
//7.2 Amortecimento
    w_meia_ind=zeros(N,N,2);      //hipermatriz com índices das freqs
                                  //da banda de meia potência
    
    w_meia_ind(:,:,1)=w_n_ind;
    w_meia_ind(:,:,2)=w_n_ind;
    
    for i=1:N                     //busca pelos pontos de meia potência
        for j=1:N
            achou_esquerda=%F;
            achou_direita=%F;
            
            H_pico=H(n,w_n_ind(j,i),j);
            H_meia_pot=abs(H_pico)/sqrt(2);
            val=[];
            
            while ~achou_esquerda
                w_meia_ind(j,i,1)=w_meia_ind(j,i,1)-1;
                
                valor_1=abs(H(n,w_meia_ind(j,i,1)-1,j));
                valor_2=abs(H(n,w_meia_ind(j,i,1),j));
                
                if (valor_1 < H_meia_pot & valor_2 > H_meia_pot) then
                    achou_esquerda=%T;
                end
            end
            
            valor_1=[];
            valor_2=[];
            
            while ~achou_direita
                w_meia_ind(j,i,2)=w_meia_ind(j,i,2)+1;
                
                valor_1=abs(H(n,w_meia_ind(j,i,2)-1,j));
                valor_2=abs(H(n,w_meia_ind(j,i,2),j));
                
                if (valor_1 > H_meia_pot & valor_2 < H_meia_pot) then
                    achou_direita=%T;
                end
            end
        end
    end
    
    for i=1:N
        for j=1:N
            etah_graf(i,j) = (w(w_meia_ind(j,i,2))^2 -...
            w(w_meia_ind(j,i,1))^2) / 2 / w_n_graf(j,i).^2;
        end
    end
    
    
    w_n_med=sum(w_n_graf,'r')/N;
    w_n_med=w_n_med';
    
    etah_val=ones(N,N);     //valores válidos para o amortecimento
                            //segundo o gráfico
    if N==6 then
      etah_val(4,1)=0;        //valores a serem desconsiderados
      etah_val(4,6)=0;
    end
    
    //cálculo da média aritmética dos valores extraídos
    etah_med=sum(etah_graf.*etah_val,'r')./sum(etah_val,'r');
    etah_med=etah_med';
    
    //matriz de autovalores experimental
    lambda_n2_med=(w_n_med.^2).*(1+%i.*etah_med);               
    
//7.3 Modos de vibrar
phi_graf=[];        //Autovetores normalizados pela massa
psi_graf=[];        //Autovalores normalizados segundo 1a linha
    
    theta=[];       //fase das ressonâcias
    
    for a=1:N       //elementos da linha 'n' de phi_graf
        for r=1:N
            A_r_jk(n,r)=imag(Hnn_nat(n,r));
            theta(n,r)=phasemag(Hnn_nat(n,r));
            
            den=%i*etah_med(r)*w_n_graf(n,r)^2/%i;
            
            phi_graf(n,r)=sqrt(A_r_jk(n,r)*den);
            phi_graf(n,r)=abs(imag(phi_graf(n,r)));
            
            phi_graf(n,r)=-theta(n,r)/abs(theta(n,r))*phi_graf(n,r);
        end
    end
    
    ind=[1:n-1,n+1:N]
    for a=ind       //restante dos elementos de phi_graf
        for r=1:N
            A_r_jk(a,r)=imag(Hnn_nat(a,r));
            theta(a,r)=phasemag(Hnn_nat(a,r));
            
            den=%i*etah_med(r)*w_n_graf(n,r)^2/%i;
            
            phi_graf(a,r)=A_r_jk(a,r)*den/phi_graf(n,r);
            phi_graf(a,r)=abs((phi_graf(a,r)));
            
            phi_graf(a,r)=-theta(a,r)/abs(theta(a,r))*phi_graf(a,r);
        end
    end
    

    
    for i=1:N
        psi_graf(:,i)=phi_graf(:,i)./phi_graf(1,i);
    end
    
    
    
//8. MODELO ESPACIAL EXPERIMENTAL
M_exp=[];               //Matriz de massa experimental
m_r_exp=[];             //Matriz de massa generalizada experimental
K_exp=[];               //Matriz de rigidez experimental
k_r_exp=[];             //Matriz de rigidez generalizada experimental
D_r_exp=[];             //Matriz de amortecimeto experimental
d_r_exp=[];             //Matriz de amortecimento generalizada
betah_exp=[];           //Constante de proporção de amortecimento 
                        //experimental
    
    m_r_exp=(inv(psi_graf)*phi_graf)^(-2); m_r_exp=real(m_r_exp);
    k_r_exp=(diag(w_n_med)^2.0)*m_r_exp;;  k_r_exp=real(k_r_exp);
    
    //eleminando elementos fora da diagnoal (muito pequenos)
    A=[];
    A=abs(m_r_exp)<(max(m_r_exp)*10^(-8));
    A=~A; A=A*1;
    m_r_exp=m_r_exp.*A;
    
    A=[];
    A=abs(k_r_exp)<(max(k_r_exp)*10^(-8));
    A=~A; A=A*1;
    k_r_exp=k_r_exp.*A;
    
    betah_exp=sum(etah_med)/N;
    d_r_exp=betah_exp*k_r_exp;
    
    M_exp=inv(psi_graf')*m_r_exp*inv(psi_graf);
    K_exp=inv(psi_graf')*k_r_exp*inv(psi_graf);
    D_exp=inv(psi_graf')*d_r_exp*inv(psi_graf);
    
    clear A
    
//9. MODELO DE RESPOSTA EXPERIMENTAL
alfa2=[];      //matriz de FRF (receptãncia)
    
    alfa2=zeros(N,N);
    
    for j=1:N
        for k=1:N
            for r=1:N
                alfa2(j,k) = alfa2(j,k)+phi_graf(j,r)*phi_graf(k,r)./ ...
                ( w_n_med(r)^2 - %z^2 + %i*etah_med(r)*w_n_med(r)^2 )
            end
        end
    end
    
    H2=zeros(N,length(w),N);
    
    for j=1:N
        for k=1:N
            H2(j,:,k)=horner(alfa2(j,k),w);
        end
    end
    
//X. EXIBIÇÃO DOS MODELOS EXPERIMENTAIS
    
    disp(' ')
    disp('[ROTA EXPERIMENTAL]');
    disp('Modelo espacial experimental');
    print(%io(2),D_exp,K_exp,M_exp);
    
    disp(' ')
    disp('Modelo modal experimental');
    print(%io(2),phi_graf,psi_graf,etah_med',w_n_med');
    
    disp(' ')
    disp('Modelo de reposta experimental (alfa_11)');
    disp(alfa2(1,1))
    
//X.1 Gráficos
     for j=1:N
         if (j>=1&j<=3) then
             scf(10)
         else
             scf(20);
         end
         subplot(211)
         plot2d('nl',w_n_graf(j,:),abs(Hnn_nat(j,:)),-7+j);
         plot2d('nl',w(w_meia_ind(j,:,1)),...
         abs(H(n,w_meia_ind(j,:,1),j)),-2);
         plot2d('nl',w(w_meia_ind(j,:,2)),...
         abs(H(n,w_meia_ind(j,:,2),j)),-2);
         
         subplot(212)
         plot2d(w_n_graf(j,:),phasemag(Hnn_nat(j,:),'m'),-N-1+j);
     end
    
    scf(30);clf();xname("FRF comparacao %d");
    subplot(211)
    plot2d('nl',w',abs([H(1,:,4)' H2(1,:,4)']),[2,3],...
            leg="H14 Original@H14 Experimental");
    subplot(212)
    plot2d(w',[phasemag([H(1,:,4); H2(1,:,4)])]',[2,3]);

//END OF FILE
