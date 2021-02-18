#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef char nev;           //részleg betűjele: a, b, c, ..., z; X-kijárat 
typedef int azon;           //részleg azonosítója: 1, 2, 3, ..., 26; 0-kijárat 
typedef int idx;            //részlegek listájában a részleg sorszáma pl: X->a->b->d->t estén d-re 3, tömb indexeléséhez hasonló



/**
 * @brief Egy részleg adatait tartalmazza
 */
typedef struct {
    azon id;                    //A részleg azonosítója (egész szám)
    int szomszedok[27];         //A részleg szomszédai: szomszedok[i]=1 ha a részleg szomszédos az i azonosítójú részleggel, 0 ha nem
    double tomeg;               //A részleg zsúfoltsága, emberek/alapterület
} reszleg;

/**
 * @brief Részlegeket tároló listaelem 
 */
typedef struct _rszlelem{
    reszleg rszl;               //Részleget tároló struktúra
    struct _rszlelem *next;     //Következő listaelem címe
}rszlelem;


/**
 * @brief Egy recept adatait tartalmazza 
 */
typedef struct{
    char nev[20+1];             //Recept neve
    int osszetevok[16];         //Recepthez szükséges összetevők: osszetevok[i]=1 ha a 2^i azonosítójú összetevő szökséges a recept elkészítéséhez, 0 ha nem
} recept;

/**
 * @brief Recepteket tároló listaelem 
 */
typedef struct _rclelem{
    recept rec;                 //Receptet tároló struktúra
    struct _rclelem *next;      //Követekező listaelem címe
} rclelem;



//------SEGÉDFÜGGVÉNYEK------//
/**
 * @brief stderr-ra kiírja a fájlkezelés/memóriafoglalás során adódó hibát
 * 
 * @param fname fájl neve
 * @param mode kezelés módja, ami során a hiba adódott: 'o' - megnyitás, 'c' - bezárás, 'm' - dinamikusan memóriafoglalás
 */
void error_message(char* fname, char mode){
    char modename[8]={0};       //"bezar" 5 karakter, "megnyit" 7 karakter
    if(mode=='m'){
        fprintf(stderr, "Sikertelen memoriafoglalas");
        return;
    }
    if(mode=='o') 
        strcpy(modename, "megnyit");
    else if(mode=='c') 
        strcpy(modename, "bezar");
    
    fprintf(stderr, "Hiba a(z) %s fajl %sasakor", fname, modename);
}

/**
 * @brief Részleg betűjelét azonosítóvá konvertálja
 * 
 * @param name részleg betűjele
 * @return azon részleg azonosítója (X->0, egyébként angol abc-beli sorszám)
 */
azon nametoid(nev name){
    if(name =='X')
        return 0;
    return name-'a'+1;
}

/**
 * @brief Kétdimenziós adjacencia mátrix koordinátáit átszámolja a tárolására használt tömb indexévé
 * 
 * @param i oszlop indexe
 * @param j sor indexe
 * @param meret mátrix mérete
 * @return int mátrix tárolására használt tömb megfelelő indexe
 */
int coord(idx i, idx j, size_t meret){
    return i*meret+j;
}

/**
 * @brief Részlegek tömbjéből a kijárat és a start elhagyásával új, permutálható tömböt hoz létre
 * 
 * @param kell részlegek sorszámainak tömbje
 * @param lk tömb hossza
 * @param perm permutálható tömb címe
 * @param start_idx  start részleg sorszáma
 * @return size_t permutálható tömb hossza
 */
size_t kellidx_to_permin(idx kell[], size_t lk, idx start_idx, idx* perm[]){
    int k=0;        //számláló
    size_t l=0;     //kimeneti tömb hossza
    for(int i=1; i<lk; i++)             //kimeneti tömb méretének meghatározása
        if(kell[i]==1 && i!=start_idx)
            l++;
    *perm = (int*)malloc(l*sizeof(int));//kimeneti tömb létrehozása
    if(*perm==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return 0;    
    }
    for(int i=1; i<lk; i++)             //és feltöltése
        if(kell[i]==1 && i!=start_idx)      //szükséges és nem a start
            (*perm)[k++]=i;
    return l;
}

/**
 * @brief Rekurzív függvény egy adott tömb elemeinek összes lehetséges sorrendjének fájlba írására
 * 
 * @param t tömb kezdőcíme
 * @param start első permutálandó elem INDEXE
 * @param stop utolsó permutálandó elem INDEXE
 * @param out kimeneti fájl címe
 */
void permutal(int t[], int start, int stop, FILE *out){
    int tmp;
    if(start==stop)
        fwrite(t, sizeof(int), stop+1, out);
    else
        for(int i=start; i<=stop; i++){
            tmp=t[start];       //csere
            t[start]=t[i];
            t[i]=tmp;
            permutal(t, start+1, stop, out); //rekurzív hívás
            tmp=t[start];       //csere vissza
            t[start]=t[i];
            t[i]=tmp;
        }
}

/**
 * @brief Kiszámolja a kapott útvonal hosszát
 * 
 * @param be útvonalat tartalmazó sorszámtömb
 * @param l tömb mérete
 * @param start_idx útvonal kezdőpontja
 * @param distm távolságokat tároló mátrix
 * @param n mátrix mérete
 * @return double útvonal hossza
 */
double uthossz(idx be[], size_t l, idx start_idx, double distm[], size_t n){
    double hossz=0;
    hossz= distm[coord(start_idx, be[0], n)];      //start és útvonal kezdőpontjának távolsága
    for(int i=1; i<l; i++)           //útvonal bejárása
        hossz += distm[coord(be[i], be[i-1], n)];    //útvonal pontjainak távolsága
    hossz += distm[coord(0, be[l-1], n)];     //útvonal végpontja és kijárat távolsága
    return hossz;
}

/**
 * @brief A kapott fájlban szereplő útvonalak közül megkeresi a legrövidebbet
 * 
 * @param utak fájl
 * @param l egy adatsor hossza a fájlban
 * @param start_idx útvonal első részlegének sorszáma
 * @param distm távolságokat tároló mátrix
 * @param n mátrix mérete
 * @return idx* minimális hosszúságú utat tartalmazó sorszámtömb
 */
idx* minutkeres(FILE* utak, size_t l, idx start_idx, double distm[], size_t n){
    double minhossz = HUGE_VAL; //minimális hosszúságú út hossza
    double behossz;             //beolvasott út hossza
    idx *minut= (idx*)malloc((l+2)*sizeof(idx));    //részlegek sorszámaiból álló minimális út
    if(minut==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return NULL;    
    }
    minut[0]=start_idx;
    minut[l+1]=0;
    idx *beUT= (idx*)malloc(l*sizeof(idx));         //részlegek sorszámaiból álló beolvasott út
    if(beUT==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return NULL;
    }

    while(fread(beUT, l, sizeof(idx), utak)!=0){
        behossz= uthossz(beUT, l, start_idx, distm, n);
        if(behossz<minhossz){           //minimumkeresés
            minhossz=behossz;
            for(int i=1, j=0; j<l; i++, j++) //útvonal bemásolása
                minut[i]=beUT[j];
        }
    }
    free(beUT);
    
    return minut;
}

//------RÉSZLEGEKEN DOLGOZÓ FÜGGVÉNYEK------//
/**
 * @brief Részlegeket tároló kétstrázsás listát hoz létre, benne a kijárattal mint segédrészleggel
 * 
 * @return rszlelem* az első strázsa címe 
 */
rszlelem* rszl_ures(void){
    rszlelem *fej= (rszlelem*)malloc(sizeof(rszlelem)); //első stázsa
    rszlelem *veg= (rszlelem*)malloc(sizeof(rszlelem)); //hátsó strázsa
    rszlelem *exit= (rszlelem*)malloc(sizeof(rszlelem)); //kijárat 
    if(fej==NULL||veg==NULL||exit==NULL){ //sikertelen foglalás
        error_message("-", 'm');
        return NULL;
        }
    fej->rszl.id=-1;        //rendezést segíti, 0 a legkisebb értékes adat
    exit->rszl.id=0; 
    exit->rszl.tomeg=0.0;   //kijáratnál nincs tömeg
    veg->rszl.id=27;        //rendezést segíti, 26 a legnagyobb értékes adat
    fej->next=exit;         //láncolás
    exit->next=veg;
    veg->next=NULL;
    return fej;
}

/**
 * @brief Új részleget hoz létre a megadott adatokból és azt az azonosítója szerint rendezve beilleszti a megadott listába
 * 
 * @param ID részleg betűjele (bemeneti fájl formátuma)
 * @param szomszedok szomszédos részlegek betűjelei szóközzel elválasztva (bemeneti fájl formátuma)
 * @param emberek részlegen tartózkodó emberek száma
 * @param terulet részleg alapterülete
 * @param fej a részlegeket tároló lista első strázsájának címe 
 */
void rszl_ujbe(nev ID, char szomszedok[51+1], int emberek, double terulet, rszlelem* fej){
    int i;
    rszlelem *rszlg=(rszlelem*)malloc(sizeof(rszlelem));    //Beszúrandó új részleg címe
    if(rszlg==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return;
    }
    //Azonosító kiszámolása
    if(ID =='X')
        rszlg->rszl.id=0;           //kijárat
    else
        rszlg->rszl.id=nametoid(ID);    //angol abc-beli sorszáma 
    //Szomszédok listájának inicializálása 
    for(i=0;i<27;i++)
        rszlg->rszl.szomszedok[i]=0;
    for(i=0;szomszedok[i]!='\n'; i++){ //és feltöltése
        int idx;
        if(szomszedok[i]==' ')      
            continue;
        else 
            idx=nametoid(szomszedok[i]);
        rszlg->rszl.szomszedok[idx]=1; 
    }
    //tomeg kiszamolasa
    rszlg->rszl.tomeg=emberek/terulet;
    //reszleg beillesztese listaba
    rszlelem *curr;
    for(curr=fej->next; curr->next->rszl.id < rszlg->rszl.id;curr=curr->next);
    rszlg->next=curr->next;
    curr->next=rszlg;
}

/**
 * @brief Egy adott betűjelű részleg címét adja meg
 * 
 * @param betu részleg betűjele (bemeneti fájl formátuma)
 * @param fej a részlegeket tároló lista első strázsájának címe
 * @return rszlelem* a megadott részleg címe, NULL ha nincs ilyen részleg
 */
rszlelem* rszl_nevtoptr(nev betu, rszlelem* fej){
    int ID= nametoid(betu);    
    for(rszlelem *curr=fej->next; curr->next; curr=curr->next)
        if(curr->rszl.id==ID)
            return curr;
    return NULL;
}

/**
 * @brief Két részleg listabeli sorszámából megadja azok távolságát (tömeg határozza meg)
 * 
 * @param r1 első részleg sorszáma
 * @param r2 második részleg sorszáma
 * @param fej a részlegeket tároló lista első strázsájának címe
 * @return double a tömeg által meghatározott távolság, ha két részleg szomszédos, egyébként végtelen 
 */
double rszl_dist(idx r1, idx r2, rszlelem* fej){
    azon id1;           //első részleg azonosítója
    double tmg1;        //első részlegben lévő tömeg
    rszlelem* p=fej->next; //pointer léptetéshez
    for(int i=0; i<r1; i++) 
        p= p->next;         //első részleg kikeresése
    id1= p->rszl.id;        //adatainak eltárolása
    tmg1= p->rszl.tomeg;
    p= fej->next;
    for(int i=0; i<r2; i++) //második részleg kikeresése
        p= p->next;
    if(p->rszl.szomszedok[id1]==0)  //nem szomszédosak
        return HUGE_VAL;
    else                            //szomszédosak
        return tmg1 +p->rszl.tomeg;
}

/**
 * @brief részlegek listájának megfelelő méretű tömböt incializál
 * 
 * @param n tömb mérete
 * @return int* részlegek feltöltetlen tömbje
 */
int* rszl_maketomb(size_t n){
    int* t =  (int*)malloc(n*sizeof(int));
    if(t==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return NULL;
    }

    t[0]=1;                 //első elem a kijárat, amire mindenképp szükség van
    for(int i=1; i<n; i++)
        t[i]=0;
    return t;
}

/**
 * @brief Az adott sorszámú részleg betűjelét adja meg
 * 
 * @param r_idx részleg sorszáma
 * @param fej a részlegeket tároló lista első strázsájának címe
 * @return nev a részleg betűjele
 */
nev rszl_idxtoname(idx r_idx, rszlelem* fej){
    azon id;                 //részleg azonosítója
    if(r_idx==0)                //kijárat
        return 'X';
    rszlelem* p=fej->next;  //pointer a léptetéshez
    for(int i=0; i<r_idx; i++)  //részleg kikeresése
        p= p->next;
    id= p->rszl.id;         
    return ('a'+id-1);
}

/**
 * @brief Felszabadítja a részlegeket tároló listát
 * 
 * @param fej a részlegeket tároló lista első strázsájának címe
 */
void rszl_freelist(rszlelem* fej){
    rszlelem* kov;      //pointer a léptetéshez
    for(kov=fej->next; fej->next; kov=kov->next){ //lista bejárása
        free(fej);                          //és felszabadítása
        fej=kov;
    }
    free(fej);
    free(kov);
}

//------RECEPTEKEN DOLGOZÓ FÜGGVÉNYEK------//
/**
 * @brief Recepteket tároló kétstrázsás üres listát hoz létre
 * 
 * @return rclelem* az első strázsa címe
 */
rclelem* rc_ures(void){
    rclelem *fej = (rclelem*)malloc(sizeof(rclelem)); //a lista első stázsája
    rclelem *veg = (rclelem*)malloc(sizeof(rclelem)); //a lista hátsó stázsája
    if(fej==NULL||veg==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return NULL;    
    }
    fej->next=veg;  //láncolás
    veg->next=NULL;
    return fej;
}

/**
 * @brief Új receptet hoz létre a megadott adatokból és azt a lista elejére szúrja
 * 
 * @param nev recept neve (bemeneti fájl formátuma)
 * @param osszetevok recepthez szükséges összetevők 1/0 sorozatként (bementei fájl formátuma)
 * @param fej a recepteket tároló lista első strázsájának címe
 */
void rc_ujbe(char nev[20+1], char osszetevok[16+1], rclelem *fej){
    rclelem *rcp=(rclelem*)malloc(sizeof(rclelem));
    if(rcp==NULL){        //sikertelen memóriafoglalás
        error_message("-", 'm');
        return;    
    }
    for(int i=0; i<20+1; i++)
        rcp->rec.nev[i]=0;
    //Név bemásolása
    strcpy(rcp->rec.nev, nev);
    //Összetevők konvertálása tömbbe
    for(int i=0; i<16; i++)
        rcp->rec.osszetevok[15-i]=osszetevok[i]-'0';
    //Lista elejáre szúrás
    rcp->next=fej->next;
    fej->next=rcp;
}

/**
 * @brief Megadott nevű recept kikeresése a listából
 * 
 * @param fej receptek listája
 * @param rcp keresett recept neve
 * @return rclelem* keresett recept listabeli címe
 */
rclelem* rc_nametoptr(rclelem* fej, char rcp[20+1]){
    rclelem* r;
    for(r=fej->next;     //szükséges recept kikeresése
        strcmp(r->rec.nev,rcp)!=0;
        r=r->next);
    return r;
}

/**
 * @brief Felszabadítja a recepteket tároló listát
 * 
 * @param fej a recepteket tároló lista első strázsájának címe
 */
void rc_freelist(rclelem* fej){
    rclelem* kov;       //pointer a léptetéshez
    for(kov=fej->next; fej->next; kov=kov->next){ //lista bejárása
        free(fej);                          //és felszabadítása
        fej=kov;
    }
    free(fej);
    free(kov);
}

//------ÖSSZETEVŐKÖN DOLGOZÓ FÜGGVÉNYEK------//
/**
 * @brief Összetevők tömbjéből törli  a recepthez nem szükségesek címeit
 * 
 * @param kell_rc bevásárolandó recept
 * @param osszetevok összetevők tömbje
 * @param ot_nevek összetevők neveit tartalmazó tömb
 * @return 0 ha sikeres, 1 ha valamelyik összetevő hiányzik
 */
int ot_removeunused(rclelem* kell_rc, rszlelem* osszetevok[], char* ot_nevek[]){
    for(int i= 0;i<16; i++){        //beolvasott összetevők címei közül
        if(kell_rc->rec.osszetevok[i]==0){   //szükségtelenek
            osszetevok[i]=NULL;             //törlése
            ot_nevek[i]=NULL;
        }
        else if(kell_rc->rec.osszetevok[i]==1 && osszetevok[i]==NULL)
            return 1;
    }
    return 0;
}

//------BEOLVASÓFÜGGVÉNYEK------//
/**
 * @brief Bevásárlóközpont adatait tartalmazó fájl egy adatsorának beolvasása a megadott változókba
 * 
 * @param in bemeneti fájl
 * @param ID részleg neve - paraméterlistán ki
 * @param szomsz szomszédok - paraméterlistán ki
 * @param fo emberek a részlegen - paraméterlistán ki
 * @param t részleg területe - paraméterlistán ki
 * @return int 1 ha sikeres, EOF ha nem
 */
int read_bevkozpont(FILE* in, nev *ID, char (*szomsz)[51+1], int *fo, double *t){
        int stat=0;
        stat= fscanf(in, "%c\n", ID);
        if(stat!=EOF){
            for(int i=0;i<51+1;i++)
                (*szomsz)[i]=0;
            fgets((*szomsz), 51, in);
            fscanf(in, "%d\n", fo);
            fscanf(in, "%lf\n", t);
        }
        return stat;
    
}

/**
 * @brief Recepteket tartalmazó fájl egy sorának beolvasása a megadott változókba
 * 
 * @param in bemeneti fájl
 * @param nev recept neve - paraméterlistán ki
 * @param otk összetevők - paraméterlistán ki
 * @return int 1 ha sikeres, EOF ha nem
 */
int read_rctek(FILE *in, char (*nev)[20+1], char (*otk)[16+1]){
    int i;                      //indexeléshez szükséges
    (*nev)[0]=fgetc(in);
    if((*nev)[0]!=EOF){
        for(i=1; (*nev)[i]=fgetc(in), (*nev)[i]!='/'; i++);
        (*nev)[i]=0;
        fgets((*otk), 17, in);
        fscanf(in, "\n");
        return 1;
    }
    return EOF;
}

/**
 * @brief Összetevőket tartalmazó fájl egy adatsorának beolvasása a megadott változókba
 * 
 * @param in bemeneti fájl neve
 * @param otidx összetevő azonosítója - paraméterlistán ki
 * @param otreszleg összetevő helye - paraméterlistán ki
 * @return int 1 ha sikeres, EOF ha nem
 */
int read_otk(FILE *in, int *otidx, nev *otreszleg){
    if(fscanf(in, "%d ", otidx)!=EOF){
        fscanf(in, "%c", otreszleg);
        return 1;
    }    
    return EOF;
}

/**
 * @brief Standard inputon érkező adatok beolvasása a megadott változókba
 *  
 * @param recept recept neve - paraméterlistán ki
 * @param start kiindulási részleg betűjele - paraméterlistán ki
 */
void read_stdin(char (*recept)[20+1], nev *start){
    int i;                      //indexeléshez kell
    printf("Adjon meg egy receptet es egy kiindulasi pontot!\n<RECEPT> <start>: ");
    for(int i=0; i<20+1; i++)
        (*recept)[i]=0;
    for(i=0; (*recept)[i]=fgetc(stdin), (*recept)[i]!=' '; i++);
    (*recept)[i]=0;
    scanf("%c", start);
}



//------A PROGRAM LOGIKAI ELEMEI------//
/**
 * @brief Bevásárlóközpont adatait tartalmazó fájl feldolgozása
 * 
 * @param fajlnev bemeneti fájl neve
 * @param num_reszlegek részlegek száma - paraméterlistán ki
 * @param fej létrehozott láncolt lista címe - paraméterlistán ki
 * @return int 0, sikertelen fájlkezelés esetén 1
 */
int process_bevkozpont(char *fajlnev, int* num_reszlegek, rszlelem* *fej){
    FILE* bevkozpont=NULL;
    nev beID;                   //részleg beolvasott betűjele
    char beSZOMSZ[51+1];        //részleg szomszédainak beolvasott felsorolása
    int beFO;                   //részlegen tartózkodó emberek beolvasott száma
    double beT;                 //részlegen beolvasott lapterülete
    *num_reszlegek=1;           //részlegek számlálója - min 1, a kijárat mint segédrészleg

    bevkozpont=fopen(fajlnev, "r");
    if(bevkozpont == NULL){
        error_message(fajlnev, 'o');
        return 1;
    }
    
    *fej=rszl_ures();

    while(read_bevkozpont(bevkozpont, &beID, &beSZOMSZ, &beFO, &beT)!=EOF){    
        rszl_ujbe(beID, beSZOMSZ, beFO, beT, *fej);
        (*num_reszlegek)++;
    }
    
    if(fclose(bevkozpont)!=0){
        error_message(fajlnev, 'c');
        return 1;
    }

    return 0;
}

/**
 * @brief Receptek adatait tartalmazó fájl feldolgozása
 * 
 * @param fajlnev bemeneti fájl neve
 * @param fej létrehozott láncolt lista címe -paraméterlistán ki
 * @return int 0, sikertelen fájlkezelés esetén 1
 */
int process_receptek(char *fajlnev, rclelem* *fej){
    FILE *rctek=NULL;
    char beRNEV[20+1];          //recept beolvasott neve
    char beOT[16+1];            //szükséges összetevők beolvasott felsorolása
    
    rctek=fopen(fajlnev, "r");
    if(rctek == NULL){
        error_message(fajlnev, 'o');
        return 1;
    }
    
    *fej=rc_ures();

    while(read_rctek(rctek, &beRNEV, &beOT)!=EOF){
        rc_ujbe(beRNEV, beOT, *fej);
    }

    if(fclose(rctek)!=0){
        error_message(fajlnev, 'c');
        return 1;
    }
    return 0;

}

/**
 * @brief Összetevőket tartalmazó fájl feldolgozása
 * 
 * @param fajlnev bemeneti fájl neve
 * @param bevkp részlegek listája
 * @param osszetevok összetevők helyét tároló tömb - paraméterlistán ki
 * @return int 0, sikertelen fájlkezelés esetén 1
 */
int process_osszetevok(char *fajlnev, rszlelem* bevkp, rszlelem *osszetevok[16]){
    FILE *otk=NULL;                  //fájlpointer
    int beOTIDX;                //összetevő beolvasott azonosítója (2^n, 0<=n<=15)
    nev beOTR;                  //összetevő helyének betűjele
    int indx;                   //beolvasott betűjelből számolt tömbbeli index

    otk=fopen(fajlnev, "r");
    if(otk == NULL){
        error_message(fajlnev, 'o');
        return 1;
    }
    
    for(int i=0;i<16;i++)
        osszetevok[i]=NULL;

    while(read_otk(otk, &beOTIDX, &beOTR)!=EOF){
        indx= round(log(beOTIDX)/log(2));
        osszetevok[indx]= rszl_nevtoptr(beOTR, bevkp);
    }
    
    if(fclose(otk)!=0){
        error_message(fajlnev, 'c');
        return 1;
    }
    return 0;
}

/**
 * @brief Standard inputon érkezett adatok feldolgozása
 * 
 * @param bevkp részlegek listája
 * @param receptek receptek listája
 * @param start kiindulási részleg azonosítója - paraméterlistán ki
 * @param recept bevásárolandó recept - paraméterlistán ki
 * @return 1 ha nem értelmezhetőek az adatok (részleg vagy recept hiányzik), különben 0
 */
int process_stdin(rszlelem* bevkp, rclelem* receptek, azon *start, char (*recept)[20+1]){
    nev beSTART;        //kiindulási részleg beolvasott betűjele
    int van=0;
    read_stdin(recept, &beSTART);
    *start=nametoid(beSTART);

    for(rszlelem* curr= bevkp->next; curr->next; curr=curr->next) //lista bejárása
        if(curr->rszl.id==*start){ //start létező részleg
            van=1;
            break;
        }
    if(!van){ //nincs
        fprintf(stderr, "A megadott reszleg nem letezik\n");
        return 1;
    }
    
    van=0;
    for(rclelem* curr= receptek->next; curr->next; curr=curr->next) //lista bejárása
        if(strcmp(*recept, curr->rec.nev)==0){ //létező recept
            van=1;
            break;
        }
    if(!van){ //nincs
        fprintf(stderr, "A megadott recept nem letezik\n");
        return 1;
    }

    return 0;
}

/**
 * @brief Létrehozza a részlegek szomszédsági mátrixát
 * 
 * @param n részlegek száma - a mátrix mérete
 * @param fej részlegek listája
 * @return double* létrehozott tömb címe
 */
double* make_adjmatrix(size_t n, rszlelem* fej){
    double *dist= (double*)malloc(n*n*sizeof(double));
    if(dist==NULL) {       //sikertelen memóriafoglalás
        error_message("-", 'm');
        return NULL;
    }
    for(int i=0; i<n; i++){     //szomszédsági mátrix feltöltése
        for(int j=0; j<n;j++)
            if(i==j)
                dist[coord(i, j, n)]=0;
            else if(j==0)
                dist[coord(i, j, n)]=rszl_dist(0, i, fej);
            else
                dist[coord(i, j, n)]=rszl_dist(i, j, fej);
    }
    return dist;
}

/**
 * @brief Floyd-algoritmust futtat a szomszédsági mátrix által megadott gráfon
 * 
 * @param dist_m mátrixot tároló tömb
 * @param n mátrix mérete 
 */
void floyd(double *dist_m, size_t n){
    for(int k=0; k<n; k++)
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            if(dist_m[coord(i, k, n)]+dist_m[coord(k, j, n)]<dist_m[coord(i, j, n)])
                dist_m[coord(i, j, n)]=dist_m[coord(i, k, n)]+dist_m[coord(k, j, n)];
}

/**
 * @brief Fájlba írj a lehetséges útvonalakat
 * 
 * @param kell_rszl részlegek tömbje
 * @param n részlegek száma (tömb mérete)
 * @param start_idx start részleg sorszáma a tömbben
 * @param l kiírt útvonalak hossza - paraméterlistán ki
 * @return int 0, sikertelen fájlkezelés esetén 1 
 */
int write_utak(int* kell_rszl, size_t n, idx start_idx, size_t *l){
    FILE *utak;     //lehetséges útvonalak eltárolására létrehozott fájl
    idx *perm_idx;  //permutáláshoz szükséges sorszámtömb címe
    
    utak=fopen("perm.bin", "wb");
    if(utak == NULL){
        error_message("perm.bin", 'o');
        return 1;
    }

    *l= kellidx_to_permin(kell_rszl, n, start_idx, &perm_idx);

    permutal(perm_idx, 0, *l-1, utak);

    if(fclose(utak)!=0){
        error_message("perm.bin", 'c');
        return 1;
    }
    utak=NULL;
    free(perm_idx);
    return 0;
}

/**
 * @brief Feldolgozza az útvonalakat tartalmazó fájlt
 * 
 * @param l egy adatsor hossza
 * @param start_idx kiindulási részleg sorszáma
 * @param distm távolságokat tároló mátrix
 * @param n mátrix mérete
 * @param minut minimális hosszúságú út tömbje
 * @return int 0, sikertelen fájlkezelés esetén 1
 */
int process_utak(size_t l, idx start_idx, double distm[], size_t n, idx *(minut[])){
    FILE* utak;

    utak=fopen("perm.bin", "rb");
    if(utak == NULL){
        error_message("perm.bin", 'o');
        return 1;
    }

    *minut= minutkeres(utak, l, start_idx, distm, n);

    if(fclose(utak)!=0){
        error_message("perm.bin", 'c');
        return 1;
    }
    utak=NULL;
    return 0;
}

/**
 * @brief felszabadítja a program által dinamikusan foglalt memóriát
 * 
 * @param rszlfej részlegek listája
 * @param rcfej receptek listája
 * @param kell_rszl részlegek tömbje
 * @param dist mátrix
 * @param minut minimális hosszúságú út tömbje
 */
void free_memory(rszlelem* rszlfej, rclelem* rcfej, int* (*kell_rszl), double* (*dist), idx* (*minut)){
    if(rszlfej)
        rszl_freelist(rszlfej);
    if(rcfej)
        rc_freelist(rcfej);
    if(*kell_rszl)
        free(*kell_rszl);
    if(*dist)
        free(*dist);
    if(*minut)
        free(*minut);
}



int main(void){
    char bevkp_nev[]="bevasarlokozpont.txt";  //részlegek adatait tartalmazó fájl neve
    char rec_nev[]="receptek.txt";              //recepteket tartalmazó fájl neve
    char otk_nev[]="osszetevok.txt";          //összetevők adatait tartalmazó fájl neve
    char *ot_nevek[16]={"liszt", "cukor", "tojas", "tej", "vaj", "sutopor", "eleszto", "vanilia", 
                        "csokolade", "mak", "lekvar", "marcipan", "mez", "fahej", "szegfuszeg", "dio"};

    int num_reszlegek;              //részlegek száma
    rszlelem* bevkp=NULL;           //részlegek listája
    rclelem* receptek=NULL;         //receptek listája
    rszlelem* osszetevok[16];       //összetevők helyét tartalmazó tömb
    azon start;                     //kiindulási részleg azonosítója
    char beRCP[20+1];               //megvásárolandó recept beolvasott neve
    double *dist=NULL;              //részlegek szomszédsági mátrixát tároló tömb (egydimenziós, kétdimenziós indexeléshez külön függvény)
    rclelem* kell_rc;               //bevásárolandó recept címe
    int *kell_rszl=NULL;            //részlegek tömbje 1-szükséges, 0-nem szükséges
    idx start_idx;                  //a start részleg helye a tömbben
    idx *minut=NULL;                //minimális hosszúságú út részlegek sorszámaként tárolva
    size_t l_kell;                  //permutáláshoz szükséges sorszámtömb mérete

    //ADATOK FELDOLGOZASA
    if(process_bevkozpont(bevkp_nev, &num_reszlegek, &bevkp)){ 
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }
    if(process_receptek(rec_nev, &receptek)){                  
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }
    if(process_osszetevok(otk_nev, bevkp, osszetevok)){
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }
    if(process_stdin(bevkp, receptek, &start, &beRCP)){
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }


    //EREDMÉNY KISZÁMÍTÁSA
    //Gráf létrehozása és legrövidebb utak
    dist= make_adjmatrix(num_reszlegek, bevkp);         //Szomszédsági mátrix létrehozása
    floyd(dist, num_reszlegek);                         //Floyd-algoritmus
    
    //Szükséges részlegek
    int j=1;        //részlegtömb indexeléséhez
    
    kell_rszl= rszl_maketomb(num_reszlegek);
    kell_rc=rc_nametoptr(receptek, beRCP);
    if(ot_removeunused(kell_rc, osszetevok, ot_nevek)){
        printf("Az adott recept osszetevoi nem talalhatok meg a bevasarlokozpontban");
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }
    for(rszlelem* curr=bevkp->next->next; 
        curr->next; 
        curr=curr->next, j++){ //részleglista és részlegtömb párhuzamos bejárása
        for(int i= 0;i<16; i++){    //összetevők helyének bejárása
            if(curr==osszetevok[i])     //aktuális részleg a szükséges címek közt
                kell_rszl[j]=1;
            if(curr->rszl.id==start){   //aktuális részleg a start
                kell_rszl[j]=1;
                start_idx=j;
            }
        }
    }

    //Legrövidebb út kiszámolása - BRUTE FORCE
    if(write_utak(kell_rszl, num_reszlegek, start_idx, &l_kell)){         
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }
    if(process_utak(l_kell, start_idx, dist, num_reszlegek, &minut)){ 
        free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
        return 1;
    }

    //Eredmény kiírása
    printf("\nBEVASARLOLISTA\n");
    for(int i=0; i<16; i++)
        if(ot_nevek[i])
            printf("- %s\n", ot_nevek[i]);
    printf("----------------\nIdealis utvonal: ");        
    for(int i=0;i<l_kell+2;i++){
        printf("%c ", rszl_idxtoname(minut[i], bevkp));
    }
    printf("\n");

    //MEMÓRIA FELSZABADÍTÁSA
    free_memory(bevkp, receptek, &kell_rszl, &dist, &minut);
    return 0;
}