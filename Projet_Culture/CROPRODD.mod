int NbParcelles =...; //nombre de parcelles
range Parcelles =1..NbParcelles;
int T =...; //horizon de planification
range Periode=1..T;
int SURF=...;
int lmax= ...; // duree au-delà de laquelle prolonger la jachère n'améliore plus le rendement'
int amax= ...; //nombre max de de semestres de cultures
{int} C[1..2] = ...; //cultures en semestre pair/impair

int s[t in Periode]=(t%2==1)?1:2;

{int} Cultures=C[1] union C[2];

 int Demande[Cultures][Periode]=...;
 
 tuple sommet{
 int l; //age de la jachere
 int a; // age de la culture
 int j; //culture ou jachere en cours 
 }
 
 {sommet} Sommets=...;
 
 tuple arc {
sommet i; //extremite initiale
sommet f; //extremite finale
int rend;
}  
 
 {arc} Arcs=...;
 {arc} InitArc=
 {<<2,0,0>,<2,0,0>,0>,
 <<2,0,0>,<2,1,1>,120>
 };
 

// Variables 
dvar boolean x[Arcs][Periode][Parcelles];
 
// Objectif
minimize sum(p in Parcelles, a in InitArc) x[a][1][p];
 
// Contraintes
subject to {
   	
    // Satisfaction de la demande
    forall (c in Cultures) {
   	    forall (t in Periode) {
            Demande[c][t] <= sum(p in Parcelles, a in Arcs : a.f.j == c) a.rend * x[a][t][p]; 
	    } 
    }

    // Contraintes de flots
    forall (p in Parcelles) {
    	sum(a in InitArc) x[a][1][p] <= 1;
   		sum(a in InitArc) x[a][1][p] == sum(a in Arcs) x[a][1][p];
	    forall (s in Sommets, t in 1..(T-1)) {
	        sum(a in Arcs : a.f == s) x[a][t][p] == sum(a in Arcs : a.i == s) x[a][t+1][p];
	    }
   }	    
}; 
