/*
	NASIOS VASILEIOS AM:3296
*/
#include<stdio.h>
#include <math.h>

float gradient[4] = {0};
float hessian[4][4] = {0};

float inv_hessian[4][4] = {  { 0.000841 ,-0.011363 , 0.043350 ,-0.041666},
							 {-0.011363 , 0.156385 ,-0.612012 , 0.607142},
							 { 0.043350 ,-0.612012 , 2.485509 ,-2.601190},
							 {-0.041666 , 0.607142 ,-2.601190 , 3.035714}	};
float p [4] = {0};  

// Pn = PNewton = -B^(-1)*gradient
void makePn(void){
	p[0] = -(inv_hessian[0][0]*gradient[0] + inv_hessian[0][1]*gradient[1] + inv_hessian[0][2]*gradient[2] + inv_hessian[0][3]*gradient[3]);
	p[1] = -(inv_hessian[1][0]*gradient[0] + inv_hessian[1][1]*gradient[1] + inv_hessian[1][2]*gradient[2] + inv_hessian[1][3]*gradient[3]);
	p[2] = -(inv_hessian[2][0]*gradient[0] + inv_hessian[2][1]*gradient[1] + inv_hessian[2][2]*gradient[2] + inv_hessian[2][3]*gradient[3]);
	p[3] = -(inv_hessian[3][0]*gradient[0] + inv_hessian[3][1]*gradient[1] + inv_hessian[3][2]*gradient[2] + inv_hessian[3][3]*gradient[3]);
}


void printP(void){
	printf("p[0] = %f , p[1] = %f , p[2] = %f , p[3] = %f\n",p[0],p[1],p[2],p[3]);
}

// compute Sum( m(t,x,y,z,w)-mt) , t from 1 to 8
float f(float x, float y, float z,float w){
	float f_total;
	float f1 = pow((x + y + z + w - 225.51) , 2);
	//printf("f1 = %f\n",f1);
	float f2 = pow((  8*x +  4*y + 2*z + w - 242.09) , 2);
	//printf("f2 = %f\n",f2);
	float f3 = pow(( 27*x +  9*y + 3*z + w - 239.23) , 2);
	//printf("f3 = %f\n",f3);
	float f4 = pow(( 64*x + 16*y + 4*z + w - 179.66) , 2);
	//printf("f4 = %f\n",f4);
	float f5 = pow((125*x + 25*y + 5*z + w - 117.60) , 2);
	//printf("f5 = %f\n",f5);
	float f6 = pow((216*x + 36*y + 6*z + w - 166.85) , 2);
	//printf("f6 = %f\n",f6);
	float f7 = pow((343*x + 49*y + 7*z + w - 181.49) , 2);
	//printf("f7 = %f\n",f7);
	float f8 = pow((512*x + 64*y + 8*z + w - 182.15) , 2);
	//printf("f8 = %f\n",f8);
	f_total = f1+f2+f3+f4+f5+f6+f7+f8;
	printf("ftotal = %f\n",f_total);
	return f_total;
}


void  makeGradient(float x, float y, float z,float w,float *gradient){

	// --- df/dx ---
	gradient[0] =    2*(x + y + z + w - 225.51)
					+ 16*(8*x +  4*y + 2*z + w - 242.09)
					+ 54*(27*x +  9*y + 3*z + w - 239.23)
					+ 128*(64*x + 16*y + 4*z + w - 179.66)
					+ 250*( 125*x + 25*y + 5*z + w - 117.60)
					+ 432*( 216*x + 36*y + 6*z + w - 166.85)
					+ 686*( 343*x + 49*y + 7*z + w - 181.49)
					+ 1024*(512*x + 64*y + 8*z + w - 182.15);

		
		

		// --- df/dy ---
	gradient[1]  =     2*(x + y + z + w - 225.51)
				      +   8*(8*x +  4*y + 2*z + w - 242.09)
	 				  + 18*(27*x +  9*y + 3*z + w - 239.23)
	 				  + 32*(64*x + 16*y + 4*z + w - 179.66)
	 				+ 50*( 125*x + 25*y + 5*z + w - 117.60) 
	                + 72*( 216*x + 36*y + 6*z + w - 166.85)
	                + 98*( 343*x + 49*y + 7*z + w - 181.49)
	                + 128*(512*x + 64*y + 8*z + w - 182.15);

	

		// --- df/dz ---
	gradient[2] =      2*(x + y + z + w - 225.51)
					+ 4*(8*x +  4*y + 2*z + w - 242.09)
					+ 6*(27*x +  9*y + 3*z + w - 239.23)
					+ 8*(64*x + 16*y + 4*z + w - 179.66)
					+ 10*(125*x + 25*y + 5*z + w - 117.60) 
					+ 12*(216*x + 36*y + 6*z + w - 166.85)
					+ 14*(343*x + 49*y + 7*z + w - 181.49)
					+ 16*(512*x + 64*y + 8*z + w - 182.15);
	

		// --- df/dw ---
	gradient[3] =     2*(x + y + z + w - 225.51)
					+ 2*(8*x +  4*y + 2*z + w - 242.09)
					+ 2*(27*x +  9*y + 3*z + w - 239.23)
					+ 2*(64*x + 16*y + 4*z + w - 179.66)
					+ 2*(125*x + 25*y + 5*z + w - 117.60) 
					+ 2*(216*x + 36*y + 6*z + w - 166.85)
					+ 2*(343*x + 49*y + 7*z + w - 181.49)
					+ 2*(512*x + 64*y + 8*z + w - 182.15);
			
		
	return;
}

void printGradient(float *gradient){

		printf("Gradient: \n df/dx = %f \n df/dy = %f \n df/dz = %f \n df/dw = %f\n",gradient[0],gradient[1],gradient[2],gradient[3]);
}

void printHessian(){
	printf("Hessian:\n");
	printf("%f  %f  %f  %f\n",hessian[0][0],hessian[0][1],hessian[0][2],hessian[0][3]);
	printf("%f  %f   %f   %f\n",hessian[1][0],hessian[1][1],hessian[1][2],hessian[1][3]);
	printf("%f   %f    %f    %f\n",hessian[2][0],hessian[2][1],hessian[2][2],hessian[2][3]);
	printf("%f    %f     %f     %f\n",hessian[3][0],hessian[3][1],hessian[3][2],hessian[3][3]);
}

void makeHessian(void){

		hessian[0][0] = 2 + 16*8 + 54*27 + 128*64 + 250*125 + 432*216 + 686*343 + 1024*512;
		hessian[0][1] = 2 + 8*8 + 18*27 + 32*64 + 50*125 + 72*216 + 98*343 + 128*512;
		hessian[0][2] = 2 + 4*8 + 6*27  + 8*64  + 10*125 + 12*216 + 14*343 + 16*512;		
		hessian[0][3] = 2 + 2*8 + 2*27  + 2*64  + 2*125  + 2*216  + 2*343  + 2*512;

		hessian[1][0] = hessian[0][1];
		hessian[1][1] = 2 + 8*4 + 18*9 + 32*16 + 50*25 + 72*36 + 98*49 + 128*64;		
		hessian[1][2] = 2 + 4*4 + 6*9  + 8*16  + 10*25 + 12*36 + 14*49 + 16*64;
		hessian[1][3] = 2 + 2*4 + 2*9  + 2*16  + 2*25  + 2*36  + 2*49  + 2*64;

		hessian[2][0] = hessian[0][2];
		hessian[2][1] = hessian[1][2];
		hessian[2][2] = 2 + 4*2 + 6*3  + 8*4  + 10*5 + 12*6 + 14*7 + 16*8;
		hessian[2][3] = 2 + 2*2 + 2*3  + 2*4  + 2*5  + 2*6  + 2*7  + 2*8;

		hessian[3][0] = hessian[0][3];
		hessian[3][1] = hessian[1][3];
		hessian[3][2] = hessian[2][3];
		hessian[3][3] = 2 + 2 + 2  + 2  + 2  + 2  + 2  + 2;
	
	return;
}

float Bisection(float a, float b){
    return (a+b)/2; 
}

// ---- Xk+1 = x + ap -----
// x[0] = x, x[1] = y, x[2] = z, x[3] =w  
void computeNext_x(float *x,float *next_x,float a){
    next_x[0] = x[0]+a*p[0];
    next_x[1] = x[1]+a*p[1];
    next_x[2] = x[2]+a*p[2];
    next_x[3] = x[3]+a*p[3];

   /* printf("next_x[0] = %f , %f , %f , %f\n",next_x[0],x[0],a,p[0]);
    printf("next_x[1] = %f , %f , %f , %f\n",next_x[1],x[1],a,p[1]);
    printf("next_x[2] = %f , %f , %f , %f\n",next_x[2],x[2],a,p[2]);
    printf("next_x[3] = %f , %f , %f , %f\n",next_x[3],x[3],a,p[3]);*/
}

float abs_ (float x){
    if(x>0){ return x;}
    return -x;
}

void lineSearch(void){
    float amax[2] = {0,2.5};        // keep ak low and ak high to reuse them
	float a[2] = {0};
    float x[4] = {-10,-10,-10,-10};	// start point is HERE!!!
    float next_x[4] = {0};
    float f1,f2;					// f1 = f(xk +ak*pk)
									// f2 = fx + c1*aj*gradientT*p
    float c1 = 0.0001;
    float c2 = 0.9;
    float aj = 0;
    float c1_a = c1*aj; 
 
    float curvature_limit;			//curvature_limit = c2*gradient(x)T*p
    float curvature;				// curvature = gradient(x+a*p)T*p
    float next_gradient[4] = {0};    
    int i =0;
    int j =0;
    
    FILE *fp;
	
	if(	(fp = fopen("Newton.txt","w")) == NULL){ return ;}

   while(i<50){ 
       
        a[0]=amax[0]; 
		a[1]=amax[1];
		
        makeGradient(x[0],x[1],x[2],x[3],gradient);
        makePn();
       
        printf("========================================================================================================\n\n");
        printf("-- i = %d\n",i);
		printf(" x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);
        printf("gradient:\n");
        printGradient(gradient);

		fprintf(fp,"%d. x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f , f = %f \n",i,x[0],x[1],x[2],x[3],f(x[0],x[1],x[2],x[3]));

        printf("Df = %f , e = 0.1\n",sqrt( gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3] ) );
        if( sqrt( gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3] ) < 0.01  ){
			printf("final x: x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);            
			printf("end\n");
            break;
        }
        printP();
        

        while(1){
            printf("\n start: j = %d \n",j);
            
            aj = Bisection(a[0],a[1]);
            printf("aj = %f\n",aj);

            computeNext_x(x,next_x, aj);
            makeGradient(next_x[0],next_x[1],next_x[2],next_x[3],next_gradient);

            printf("next x: next_x[0] = %f , next_x[1] = %f , next_x[2] = %f , next_x[3] = %f\n",next_x[0],next_x[1],next_x[2],next_x[3]);
            printf("next_gradient:\n");
            printGradient(next_gradient);

			// f1 = f(xk +ak*pk)
            f1 = f(next_x[0],next_x[1],next_x[2],next_x[3]);
			// f2 = fx + c1*aj*gradientT*p                  
            f2 = f(x[0],x[1],x[2],x[3]) + c1*aj*(gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3]);        
            printf("f1 = %f\n",f1);
            printf("f2 = %f , c1*aj*... = %f , c1 = %f , aj = %f , ... = %f\n",f2,c1*aj*(gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3]),c1,aj,gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3]);   
    
            if(f1>=f2 && abs_(f1-f2)>0.01){  //--- Armijo
                printf("i am in if\n");
				printf("f1-f2 = %f\n",abs_(f2-f1));
                a[1] = aj;
                printf("a[0] = %f , a[1] = %f\n",a[0],a[1]);
            }else{
                printf("i am in else\n");
                // curvature = gradient(x+a*p)T*p
                curvature = next_gradient[0]*p[0] + next_gradient[1]*p[1] + next_gradient[2]*p[2] + next_gradient[3]*p[3];                
				
				//curvature_limit = c2*gradient(x)T*p
				curvature_limit = c2*( gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3] ); 
                printf("curvature = %f , curvature_limit = %f \n",curvature,curvature_limit);


                
                if( (curvature>=curvature_limit && curvature<=-curvature_limit) || curvature==curvature_limit){       
                    printf("wolfe conditions ok\n");
					printf("curvature = %f , curvature_limit = %f ,  curvature_limit = %f \n",curvature,curvature_limit,-curvature_limit);
                    x[0] = next_x[0]; x[1] = next_x[1]; x[2] = next_x[2]; x[3] = next_x[3];
                    printf("next x: x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);
                    break;
                }
                if(curvature > -curvature_limit ){
                    printf("curvature > -curvature_limit = %f>%f => up\n",curvature,-curvature_limit);
					a[1] = aj;                    
                    printf("a[0] = %f , a[1] = %f\n",a[0],a[1]);
                }
                if(curvature < curvature_limit){
                    printf("curvature < curvature_limit = %f<%f =>down\n",curvature,curvature_limit);
                    a[0] = aj;
                    printf("a[0] = %f , a[1] = %f\n",a[0],a[1]);
                }              
            }
            j++; 
        }
        i++;  
        j=0;      
    }
	fclose(fp);
}

// compute B*gradient
void multiplyBk_gk(float B_g[]){
	B_g[0] = (hessian[0][0]*gradient[0] + hessian[0][1]*gradient[1] + hessian[0][2]*gradient[2] + hessian[0][3]*gradient[3]);
	B_g[1] = (hessian[1][0]*gradient[0] + hessian[1][1]*gradient[1] + hessian[1][2]*gradient[2] + hessian[1][3]*gradient[3]);
	B_g[2] = (hessian[2][0]*gradient[0] + hessian[2][1]*gradient[1] + hessian[2][2]*gradient[2] + hessian[2][3]*gradient[3]);
	B_g[3] = (hessian[3][0]*gradient[0] + hessian[3][1]*gradient[1] + hessian[3][2]*gradient[2] + hessian[3][3]*gradient[3]);	
}

/*   getDirection is used by dogleg  */
void getDirection(float x[],float radious){

	float gT_B_g;
	float B_g[4] ={0};
	float division;					// division = - gT*g/gT*B*g
	float pU[4] = {0};  		    // pU = division * gardient
	float add[4]={0};				// add = pU + (t-1)(p-pU)
	float p_norm,gradient_norm;
	float gT_g;
	float pU_norm;

	makeGradient(x[0],x[1],x[2],x[3],gradient);
	//printGradient(gradient);
	makeHessian();
	makePn();
	//printP();

	p_norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
 	printf("p_norm = %f , radious = %f\n",p_norm,radious);

	if( p_norm <= radious){  // norm_p<=radious for dogleg
		printf("dogleg: p[0] = %f , p[1] = %f , p[2] = %f , p[3] = %f \n",p[0],p[1],p[2],p[3]);
		return;
	}

	gradient_norm = sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3]);
	gT_g = gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3];		
	multiplyBk_gk(B_g);
	printf("B_g[0] = %f , B_g[1] = %f , B_g[2] = %f , B_g[3] = %f \n",B_g[0],B_g[1],B_g[2],B_g[3]);
	gT_B_g = gradient[0]*B_g[0] + gradient[1]*B_g[1] + gradient[2]*B_g[2] + gradient[3]*B_g[3];
	division = - gT_g/gT_B_g;

	printf("division = %f , gT_g = %f , gT_B_g = %f \n",division,gT_g,gT_B_g);
	

	pU[0] = division*gradient[0];
	pU[1] = division*gradient[1];
	pU[2] = division*gradient[2];
	pU[3] = division*gradient[3];
 		
	pU_norm = sqrt( pU[0]*pU[0] + pU[1]*pU[1] + pU[2]*pU[2] + pU[3]*pU[3]);	
	printf("pU: pU[0] = %f , pU[1] = %f , pU[2] = %f , pU[3] = %f \n",pU[0],pU[1],pU[2],pU[3]);	
	printf("pU_norm = %f , radious = %f\n",pU_norm,radious);

	if( pU_norm >= radious ){  //  PU_norn >= radious  for cauchy
		p[0] = -(radious/gradient_norm)*gradient[0];
		p[1] = -(radious/gradient_norm)*gradient[1];
		p[2] = -(radious/gradient_norm)*gradient[2];
		p[3] = -(radious/gradient_norm)*gradient[3];
		printf("cauchy point: p[0] = %f , p[1] = %f , p[2] = %f , p[3] = %f \n",p[0],p[1],p[2],p[3]);
		return;
	}

	//  pU_norm < radious  solve  s = [p U + (τ − 1)(p B − p U )]^2 - radious^2   
	//  t is in [1,2]
	float t[2] = {1,2};
	float t2;
	float s;
	float r = radious*radious;
	
	while(1){
		t2 = (t[1] + t[0])/2;//+t[0];
		add[0] = pU[0] + (t2-1)*(p[0]-pU[0]);
		add[1] = pU[1] + (t2-1)*(p[1]-pU[1]);
		add[2] = pU[2] + (t2-1)*(p[2]-pU[2]);
		add[3] = pU[3] + (t2-1)*(p[3]-pU[3]);
		s = add[0]*add[0] + add[1]*add[1] + add[2]*add[2] + add[3]*add[3] - r;
		printf("s = %f , t2 = %f , t[0] = %f, t[1] = %f\n",s,t2,t[0],t[1]);
		//printf("abs(s) = %f\n",abs_(s));

		if(s==0 || abs_(s)<0.0001 || t[0]==t2 || t[1]==t2){ printf("i am in break\n"); break; }
		if(s>0){t[1]=t2;}
		if(s<0){t[0]=t2;}
	}

	p[0] = add[0];
	p[1] = add[1];
	p[2] = add[2];
	p[3] = add[3];
	printf("cauchy2 point: p[0] = %f , p[1] = %f , p[2] = %f , p[3] = %f \n",p[0],p[1],p[2],p[3]);				
}

// squere model:  m(p) = f + gT*p + 1/2*p_T*B*p
// input x and compute f , gradient
float m(float x[]){
	float B_p[4] ={0};
	float f1 = f(x[0],x[1],x[2],x[3]);	
	float gT_p = gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3];
	

	B_p[0] = (hessian[0][0]*p[0] + hessian[0][1]*p[1] + hessian[0][2]*p[2] + hessian[0][3]*p[3]);
	B_p[1] = (hessian[1][0]*p[0] + hessian[1][1]*p[1] + hessian[1][2]*p[2] + hessian[1][3]*p[3]);
	B_p[2] = (hessian[2][0]*p[0] + hessian[2][1]*p[1] + hessian[2][2]*p[2] + hessian[2][3]*p[3]);
	B_p[3] = (hessian[3][0]*p[0] + hessian[3][1]*p[1] + hessian[3][2]*p[2] + hessian[3][3]*p[3]);
	
	float pT_B_p = p[0]*B_p[0] + p[1]*B_p[1] + p[2]*B_p[2] + p[3]*B_p[3];
	printf("m=%f , f1=%f , gT_p=%f , 0.5*pT_B_p=%f\n",f1 + gT_p + 0.5*pT_B_p,f1,gT_p,0.5*pT_B_p);

	return (f1 + gT_p + 0.5*pT_B_p);
}

void dogleg(){
	float x[4] = {1,1,1,1};		// srart point is HERE !!!
	float radious_max = 10;
	float radious = 5;
	float m_0;
	float df,dm,r,p_norm;   	// df = fk - fk+p  , dm = m0 - mx , r = dx/dm
	

	int i=0;
	FILE *fp;	
	if(	(fp = fopen("Dogleg.txt","w")) == NULL){ return ;}

	while(i<30){

		printf("\n -------------------- i=%d ------------------\n",i);
		printf("x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);

		makeGradient(x[0],x[1],x[2],x[3],gradient);
		makePn();

		fprintf(fp,"%d. x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f , f = %f \n",i,x[0],x[1],x[2],x[3],f(x[0],x[1],x[2],x[3]));
		printf("Df = %f , e = 0.1\n",sqrt( gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3] ) );

		// if ||gradient|| < 0.1 stop
        if( sqrt( gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3] ) < 0.1  ){
			printf("final x: x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);            
			printf("end\n");
            break;
        }

		getDirection(x,radious);
		//printf("get direction: p[0] = %f , p[1] = %f , p[2] = %f, p[3] = %f \n",p[0],p[1],p[2],p[3]);
		
		// compute r = df/dm
		p_norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
		m_0 = f(x[0],x[1],x[2],x[3]);
		df = f(x[0],x[1],x[2],x[3]) - f(x[0]+p[0],x[1]+p[1],x[2]+p[2],x[3]+p[3]);
		dm = m_0 - m(x);	
		r = df/dm;
		printf("r=%f , df=%f , dm=%f\n",r,df,dm);

		/*if(r<0 || r>1){
			continue;
		}*/

		if(r<0.25){
			radious=0.25*radious;
		}
		if(r>0.75 /*&& p_norm==radious*/){
			printf("2*radious=%f , radious_max=%f\n",2*radious,radious_max);
			if((2*radious)>=radious_max){
				radious = radious_max;
			}else{ 
				printf("hi2\n");
			} 
			printf("radious=%f\n",radious);
		}
		if(r>0.2){
			printf("r>0.2 change x\n\n\n");
			x[0] = x[0] + p[0];
			x[1] = x[1] + p[1];
			x[2] = x[2] + p[2];
			x[3] = x[3] + p[3];
		}	
		i++;
		
	}
	fclose(fp);
}

// H = B^-1  used by BFGS 
float H[4][4] = {{0.5,0,0,0},{0,0.5,0,0},{0,0,0.5,0},{0,0,0,0.5}};

// compute p for BFGS
void makePBFGS(void){
	p[0] = -(H[0][0]*gradient[0] + H[0][1]*gradient[1] + H[0][2]*gradient[2] + H[0][3]*gradient[3]);
	p[1] = -(H[1][0]*gradient[0] + H[1][1]*gradient[1] + H[1][2]*gradient[2] + H[1][3]*gradient[3]);
	p[2] = -(H[2][0]*gradient[0] + H[2][1]*gradient[1] + H[2][2]*gradient[2] + H[2][3]*gradient[3]);
	p[3] = -(H[3][0]*gradient[0] + H[3][1]*gradient[1] + H[3][2]*gradient[2] + H[3][3]*gradient[3]);
}


void printH(){
	int i,j;
	printf("H:\n");
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			printf("%f  ",H[i][j]);	
		}
		printf("\n");
	}
	printf("\n");
}

void BFGS(){

	float x[4] = {1,1,1,1};			// start point is HERE!!!
	float next_x[4] = {0};
	float next_gradient[4] = {0};
	float s[4] = {0};				// s = Xk+1 - Xk
	float y[4] = {0};				// y = gradient(Xk+1) - gradient(Xk)
		float r;					// r = i/(yT*s) 
	float r_s[4] = {0};				// r_s = r*s
	float r_y[4] = {0};				// r_y = r*y


	float s1[4][4] = {0}; 			// s1 = (I − r*s*yT)H 
	float s2[4][4] = {0}; 			// s2 = s1*(I − ρ*y*sT)

	int i,j,k;

//----- the same with line search
    float amax[2] = {0,0.001};
	float a[2] = {0};
    float f1,f2;
    float c1 = 0.0001;
    float c2 = 0.9;
    float aj = 0;
    float c1_a = c1*aj; 
 
    float curvature_limit;
    float curvature;  
	int counter=0;
	int counter2=0;
    
    FILE *fp;
	
	if(	(fp = fopen("BFGS.txt","w")) == NULL){ return ;}
		
	while(counter<20){
        printf("--------------------------- counter = %d  ----------------------------------\n",counter);
					
        a[0]=amax[0]; a[1]=amax[1];
		//gradient[0]=0; gradient[1]=0; gradient[2]=0; gradient[3]=0;

        makeGradient(x[0],x[1],x[2],x[3],gradient);
        makePBFGS();
 
       
		printf(" x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);
        printf("gradient:\n");
        printGradient(gradient);

		fprintf(fp,"%d. x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f , f = %f \n",counter,x[0],x[1],x[2],x[3],f(x[0],x[1],x[2],x[3]));

        printf("Df = %f , e = 0.1\n",sqrt( gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3] ) );
        if( sqrt( gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2] + gradient[3]*gradient[3] ) < 0.01  ){
			printf("final x: x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);            
			printf("end\n");
            break;
        }
        
        counter2=0;
		while(counter2<10){
            printf("\n start: counter2 = %d \n",counter2);
            
            aj = Bisection(a[0],a[1]);
            printf("aj = %f\n",aj);

            computeNext_x(x,next_x, aj);
            makeGradient(next_x[0],next_x[1],next_x[2],next_x[3],next_gradient);

            printf("next x: next_x[0] = %f , next_x[1] = %f , next_x[2] = %f , next_x[3] = %f\n",next_x[0],next_x[1],next_x[2],next_x[3]);
            printf("next_gradient:\n");
            printGradient(next_gradient);

			// f1 = f(x+a*p)
            f1 = f(next_x[0],next_x[1],next_x[2],next_x[3]);                
			// f2 = f(x) + c1*a*gardientT*p            
			f2 = f(x[0],x[1],x[2],x[3]) + c1*aj*(gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3]);        
            printf("f1 = %f\n",f1);
            printf("f2 = %f , c1*aj*... = %f , c1 = %f , aj = %f , ... = %f\n",f2,c1*aj*(gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3]),c1,aj,gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3]);   
    
            if(f1>=f2 && abs_(f1-f2)>0.01){  //f1>=f2   Armijo: f1 must be < f2
                printf("i am in if\n");
				printf("f1-f2 = %f\n",abs_(f2-f1));
                a[1] = aj;
                printf("a[0] = %f , a[1] = %f\n",a[0],a[1]);
            }else{
                printf("i am in else\n");
                
				// curvature = gradient(x+a*p)T*p
                curvature = next_gradient[0]*p[0] + next_gradient[1]*p[1] + next_gradient[2]*p[2] + next_gradient[3]*p[3];
				// corvature_limit = c2*gradient(x)T*p                 
				curvature_limit = c2*( gradient[0]*p[0] + gradient[1]*p[1] + gradient[2]*p[2] + gradient[3]*p[3] ); 
                printf("curvature = %f , curvature_limit = %f \n",curvature,curvature_limit);

                
                if( (curvature>=curvature_limit && curvature<=-curvature_limit) || curvature==curvature_limit){       // strong conditions
                    printf("wolfe conditions ok\n");
					printf("curvature = %f , curvature_limit = %f ,  curvature_limit = %f \n",curvature,curvature_limit,-curvature_limit);
                    printf("next x: x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);
                    break;
                }
                if(curvature > -curvature_limit ){
                    printf("curvature > -curvature_limit = %f>%f => up\n",curvature,-curvature_limit);
					a[1] = aj;                    
                    printf("a[0] = %f , a[1] = %f\n",a[0],a[1]);
                }
                if(curvature < curvature_limit){
                    printf("curvature < curvature_limit = %f<%f =>down\n",curvature,curvature_limit);
                    a[0] = aj;
                    printf("a[0] = %f , a[1] = %f\n",a[0],a[1]);
                }              
            } 
			counter2++;
        }
		


//--------------------------------------------	 
		printf("\nstart making H:\n");

		s[0] = next_x[0] - x[0];
		s[1] = next_x[1] - x[1];
		s[2] = next_x[2] - x[2];
		s[3] = next_x[3] - x[3];
		printf("x: x[0] = %f , x[1] = %f , x[2] = %f, x[3] = %f \n",x[0],x[1],x[2],x[3]);
		printf("next_x: next_x[0] = %f , next_x[1] = %f , next_x[2] = %f, next_x[3] = %f \n\n",next_x[0],next_x[1],next_x[2],next_x[3]);
		printf("s[0] = %f , s[1] = %f , s[2] = %f ,s[3] = %f\n",s[0],s[1],s[2],s[3]);

		// next_x is valid so x = Xk+1
		x[0] = next_x[0]; x[1] = next_x[1]; x[2] = next_x[2]; x[3] = next_x[3];


		y[0] = next_gradient[0]-gradient[0];
		y[1] = next_gradient[1]-gradient[1];
		y[2] = next_gradient[2]-gradient[2];
		y[3] = next_gradient[3]-gradient[3];
		printf("y[0] = %f , y[1] = %f , y[2] = %f ,y[3] = %f\n",y[0],y[1],y[2],y[3]);
		
		r = 1/(y[0]*s[0] + y[1]*s[1] + y[2]*s[2] + y[3]*s[3]);
		printf("r = %f\n",r);
		
		r_s[0] = r*s[0];
		r_s[1] = r*s[1];
		r_s[2] = r*s[2];
		r_s[3] = r*s[3];
		printf("r_s[0] = %f , r_s[1] = %f , r_s[2] = %f ,r_s[3] = %f\n",r_s[0],r_s[1],r_s[2],r_s[3]);

		r_y[0] = r*y[0];
		r_y[1] = r*y[1];
		r_y[2] = r*y[2];
		r_y[3] = r*y[3];
		printf("r_y[0] = %f , r_y[1] = %f , r_y[2] = %f ,r_y[3] = %f\n",r_y[0],r_y[1],r_y[2],r_y[3]);
		
		float r_s_yT[4][4] = {  { 1-r_s[0]*y[0] , -r_s[0]*y[1]   , -r_s[0]*y[2]  , -r_s[0]*y[3]} ,
								  { -r_s[1]*y[0]  , 1-r_s[1]*y[1]  , -r_s[1]*y[2]  , -r_s[1]*y[3]} ,
								  { -r_s[2]*y[0]  , -r_s[2]*y[1]   , 1-r_s[2]*y[2] , -r_s[2]*y[3]} ,
								  { -r_s[3]*y[0]  , -r_s[3]*y[1]   , -r_s[3]*y[2]  , 1-r_s[3]*y[3]}    }; 
	
		float r_y_sT[4][4] = {    { 1-r_y[0]*s[0] , -r_y[0]*s[1]  , -r_y[0]*s[2]  , -r_y[0]*s[3] } ,
									{ -r_y[1]*s[0]  , 1-r_y[1]*s[1] , -r_y[1]*s[2]  , -r_y[1]*s[3] } ,
									{ -r_y[2]*s[0]  , -r_y[2]*s[1]  , 1-r_y[2]*s[2] , -r_y[2]*s[3] } ,
									{ -r_y[3]*s[0]  , -r_y[3]*s[1]  , -r_y[3]*s[2]  , 1-r_y[3]*s[3] }   };
		
		float r_s_sT[4][4] = {	{ r_s[0]*s[0] , r_s[0]*s[1] , r_s[0]*s[2] , r_s[0]*s[3]} ,
								{ r_s[1]*s[0] , r_s[1]*s[1] , r_s[1]*s[2] , r_s[1]*s[3]} ,
								{ r_s[2]*s[0] , r_s[2]*s[1] , r_s[2]*s[2] , r_s[2]*s[3]} ,
								{ r_s[3]*s[0] , r_s[3]*s[1] , r_s[3]*s[2] , r_s[3]*s[3]}   };

	



		// compute s1
		for(i=0; i<4; i++){
			for(j=0; j<4; j++){
				for(k=0; k<4; k++){
					s1[i][j] += r_s_yT[i][k]*H[k][j];
				}
			}
		}

		// compute s2
		for(i=0; i<4; i++){
			for(j=0; j<4; j++){
				for(k=0; k<4; k++){
					s2[i][j] += s1[i][k]*r_y_sT[k][j];
				}
			}
		}



		//	compute new H	
		for(i=0; i<4; i++){
			for(j=0; j<4; j++){
				H[i][j] = s2[i][j] + r_s_sT[i][j];
			}
		}

// ---  now print all matrices ---
	printf("r_s_yT:\n");
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			printf("%f  ",r_s_yT[i][j]);	
		}
		printf("\n");
	}
	printf("\n");

	printf("r_y_sT:\n");
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			printf("%f  ",r_y_sT[i][j]);	
		}
		printf("\n");
	}
	printf("\n");

	printf("r_s_sT:\n");
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			printf("%f  ",r_s_sT[i][j]);	
		}
		printf("\n");
	}
	printf("\n");

	printf("s1:\n");
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			printf("%f  ",s1[i][j]);	
		}
		printf("\n");
	}
	printf("\n");

	printf("s2:\n");
	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			printf("%f  ",s2[i][j]);	
		}
		printf("\n");
	}
	printf("\n");


		printH();
		counter++;
	}
	fclose(fp);
}




int main(void){
    lineSearch();	
	dogleg();	
	BFGS();
}
