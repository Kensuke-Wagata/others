
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mt19937.h"

long Get_Distance(double, double, double, double);
double Get_Circuit_Cost(long**, int*, double*, double*);
void GNU_Plot(int, int*, double*, double*);
double Random_Search(long**, double*, double*);
double Hill_Random_Exchange(long**, double*, double*);
double Hill_2_Opt(long**, double*, double*);
double Simulated_Annealing(long**, double*, double*);
void Get_Route_Array(int*);
double getrusage_sec();

static int  PLOT_INTERVAL = 300000000;//プロット時の回数の間隔
static int  INIT_T;
static int INIT_R;
static double  ALFA;

static int n;
static int times;
static char *data_file;
static char *solve_method;
static char *vicinity;
static int hill_MODE;
static int MODE;
static int try_num = 0;
static double execution_time = 0;
static int SA_ans = 0;
static int set_seed = 0;

int main( int argc, char *argv[] ){
    int i, j, input;
    double min, result;
    FILE *fp, *para_fp;
    char temp[100];
    double parameter;
    
    //コマンドライン入力確認
    if( argc!=3 ){
        printf( "Usa1ge: sample <input_filename> <execution_mode>\n" );
        exit( 1 );
    }
    if(( fp=fopen( argv[1], "r" ))==NULL ){
        printf( "file open error!\n" );
        exit( 1 );
    }
    if(strcmp(argv[2], "random") != 0 &&
       strcmp(argv[2], "hill-climbing") != 0 &&
       strcmp(argv[2], "simulated-annealing") != 0) {
        printf( "execution mode error!\n" );
        printf( "random | hill-climbing | simulated-annealing\n" );
        exit( 1 );
    }
    
    //コマンドライン入力に対しての手法、近傍の設定
    data_file = argv[1];
    solve_method = argv[2];
    if(strcmp(argv[2], "random") == 0 ) MODE = 1;
    if(strcmp(argv[2], "hill-climbing") == 0) MODE = 2;
    if(strcmp(argv[2], "simulated-annealing") == 0) MODE = 3;
    if(MODE == 2){
        printf("山登り法の近傍（1->ランダム交換、2->2-opt）：");
        scanf("%d", &hill_MODE);
        if(hill_MODE == 1) vicinity = "ランダム選択";
        else vicinity = "2-opt";
    }
    
    //パラメータファイル読み込み
    if(MODE == 3) {
        vicinity = "2-opt";
        if(( para_fp=fopen( "parameter_file.txt", "r" ))==NULL ){
            printf( "parameter_file open error!\n" );
            exit( 1 );
        }
        
        do{
            fscanf(para_fp, "%s %lf", temp, &parameter);
            if(strcmp(temp, "INIT_T") == 0) INIT_T = (int)parameter;
            if(strcmp(temp, "INIT_R") == 0) INIT_R = (int)parameter;
            if(strcmp(temp, "ALFA") == 0) ALFA = parameter;
        }while(strcmp(temp, "EOF") != 0);
        fclose(para_fp);
    }
    
    //データファイルで数字データ部分まで読み込み
    do {
        fscanf( fp, "%s", temp );
        if( strcmp( "DIMENSION", temp )==0 ){
            fscanf( fp, "%s", temp );
            break;
        }
        if( strcmp( "DIMENSION:", temp )==0 ) break;
    } while( 1 );
    
    fscanf( fp, "%d", &n);
    fclose( fp );
    
    //---------------------------------------------------------------
    
    double *x = (double*)malloc(sizeof(double) * n);
    double *y = (double*)malloc(sizeof(double) * n);
    long **distance = (long**)malloc(sizeof(long) * n);
    int trash;
    
    for(i = 0; i < n; i++)
        distance[i] = (long*)malloc(sizeof(long) * n);
    
    //ファイルの数字データ部分まで読み込み
    if(( fp=fopen( argv[1], "r" ))==NULL ){
        printf( "file open error!\n" );
        exit( 1 );
    }
    while(1){
        fscanf(fp, "%s", temp);
        if(strcmp("NODE_COORD_SECTION", temp) == 0)
            break;
    }
    // 数字読み込み
    for(i = 0; i < n; i++){
        fscanf(fp, "%d %lf %lf", &trash, &x[i], &y[i]);
    }
    //各都市間の距離格納
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            distance[i][j] = Get_Distance(x[i], y[i], x[j], y[j]);
        }
    }
    
    //結果書き込み用ファイル準備
    FILE *data_fp;
    char *filename = "SA_usa500000.txt";
    if( (data_fp = fopen(filename, "w") ) == NULL){
        printf("File Open Error: %s\n", filename);
        exit(1);
    }
    
    times = 1;
    for(i = 0; i < 1000; i++){
        if(MODE == 1){//ランダム探索
            min = Random_Search(distance, x, y);
            fprintf(data_fp, "%d %lf\n", i+1, min);
            printf("%d %lf\n", i+1, min);
            times++;
        }
        if(MODE == 2){//山登り法
            if(hill_MODE == 1) {//近傍：ランダム交換
                min = Hill_Random_Exchange(distance, x, y);
                fprintf(data_fp, "%d %lf\n", i+1, min);
                printf("%d %lf\n", i+1, min);
                times++;
            }
            if(hill_MODE == 2) {//近傍：2-opt
                min = Hill_2_Opt(distance, x, y);
                fprintf(data_fp, "%d %lf\n", i+1, min);
                printf("%d %lf\n", i+1, min);
                times++;
            }
        }
        if(MODE == 3){//シミュレーティッド・アニーリング法
            min = Simulated_Annealing(distance, x, y);
            fprintf(data_fp, "%d %lf\n", i+1, min);
            printf("%d %lf\n", i+1, min);
            times++;
        }
        try_num = 0;
    }
    
    fclose(data_fp);
    fclose( fp );
}


//ランダム探索
double Random_Search(long **distance, double x[], double y[]){
    int i;
    double min, result;
    int route[n];
    for(i = 0; i < times; i++){
        Get_Route_Array(route);
        if(i == 0) min = Get_Circuit_Cost(distance, route, x, y);
        else {
            result = Get_Circuit_Cost(distance, route,  x, y);
            if(result < min) min = result;
        }
    }
    return min;
}


//山登り法
double Hill_Random_Exchange(long **distance, double x[], double y[]){
    int route[n];
    int init, seed, i, j, tmp, num;
    int point[2];
    int best_point[2];
    double random, result, min;
    bool flag = true;
    
    //初期ルートを取得、コストを計算
    Get_Route_Array(route);
    init = Get_Circuit_Cost(distance, route, x, y);
    min = init;
    
    num = n * (n-1) / 2;
    
    while(flag){
        flag = false;
        init = min;
        sgenrand(times);
        for(i = 0; i < num; i++){
            //２つの都市を選ぶ
            random = genrand()*1000;
            point[0] = (int)random;
            point[0] %= n;
            random = genrand()*1000;
            point[1] = (int)random;
            point[1] %= n;
            
            //ルートをswap
            tmp = route[point[0]];
            route[point[0]] = route[point[1]];
            route[point[1]] = tmp;
            
            //ルートのコスト計算
            result = Get_Circuit_Cost(distance, route, x, y);
            
            //routeを戻す
            tmp = route[point[0]];
            route[point[0]] = route[point[1]];
            route[point[1]] = tmp;
            
            if(result < init){
                flag = true;
                if(result < min) {//探索中の最小暫定解時の2つの都市データを記憶
                    min = result;
                    best_point[0] = point[0];
                    best_point[1] = point[1];
                }
            }
        }
        //route更新
        tmp = route[best_point[0]];
        route[best_point[0]] = route[best_point[1]];
        route[best_point[1]] = tmp;
    }
    return min;
}


//山登り法2-opt
double Hill_2_Opt(long **distance, double x[], double y[]){
    int route[n];
    int seed, i, j, p, tmp;
    int num, element, swap_num;
    int point[2];
    int best_point[2];
    double random, result, min, init;
    bool flag = true;
    int count;
    
    //初期ルートを取得、コストを計算
    Get_Route_Array(route);
    init = Get_Circuit_Cost(distance, route, x, y);
    min = init;
    
    num = n * (n-1) / 2;
    
    while(flag){
        flag = false;
        init = min;
        sgenrand(times);
        for(i = 0; i < num; i++){
            //２つの都市を選ぶ
            random = genrand()*1000;
            point[0] = (int)random;
            point[0] %= n;
            random = genrand()*1000;
            point[1] = (int)random;
            point[1] %= n;
            //point[0] < point[1]なるよう入れ替え
            if(point[0] > point[1]){
                tmp = point[0];
                point[0] = point[1];
                point[1] = tmp;
            }
            
            //スワップ回数の計算
            element = point[1] - point[0] + 1;
            if(element % 2 != 0) swap_num = (element-1)/2;
            if(element % 2 == 0) swap_num = element/2;
            if(element == 0) continue;
            
            //routeをスワップ
            count = 0;
            for(p = point[0]; p < point[0] + swap_num; p++){
                tmp = route[p];
                route[p] = route[ point[1]-count ];
                route[ point[1]-count ] = tmp;
                count++;
            }
            
            //ルートのコストを計算
            result = Get_Circuit_Cost(distance, route, x, y);
            
            //route戻す
            count = 0;
            for(p = point[0]; p < point[0] + swap_num; p++){
                tmp = route[p];
                route[p] = route[ point[1]-count ];
                route[ point[1]-count ] = tmp;
                count++;
            }
            
            if(result < init){
                flag = true;
                if(result < min) {//探索中の最小暫定解時の2つの都市データを記憶
                    min = result;
                    best_point[0] = point[0];
                    best_point[1] = point[1];
                }
            }
        }
        //route更新
        element = best_point[1] - best_point[0] + 1;
        if(element % 2 != 0) swap_num = (element-1)/2;
        if(element % 2 == 0) swap_num = element/2;
        count = 0;
        for(p = best_point[0]; p < best_point[0] + swap_num; p++){
            tmp = route[p];
            route[p] = route[ best_point[1]-count ];
            route[ best_point[1]-count ] = tmp;
            count++;
        }
    }
    return min;
}


//シミュレーティッド・アニーリング法
double Simulated_Annealing(long **distance, double x[], double y[]){
    int route[n];
    int seed, i, j, p, tmp;
    int num, element, swap_num, D;
    int point[2];
    double random, init, result, min;
    int count;
    double T = INIT_T;
    double R = INIT_R;
    bool flag = true;
    int n_t = 0;
    double t1, t2;
    
    //初期ルート、コストを計算
    Get_Route_Array(route);
    init = Get_Circuit_Cost(distance, route, x, y);
    min = init;
    
    num = n * (n-1) / 2;
    
    //計算時間計測スタート
    t1 = 0; t2 = 0; n_t = 0; execution_time = 0;
    t1 = getrusage_sec();
    n_t *= i;
    
    while(T > INIT_T/20 && flag){
        T = INIT_T;
        R = INIT_R;
        
        for(i = 0; i < R; i++){
            flag = false;
            sgenrand(times + try_num);
            
            //２つの都市を選ぶ
            random = genrand()*1000;
            point[0] = (int)random;
            point[0] %= n;
            random = genrand()*1000;
            point[1] = (int)random;
            point[1] %= n;
            //point[0] < point[1]になるよう入れ替え
            if(point[0] > point[1]){
                tmp = point[0];
                point[0] = point[1];
                point[1] = tmp;
            }
            
            //スワップ回数を計算
            element = point[1] - point[0] + 1;
            if(element % 2 != 0) swap_num = (element-1)/2;
            if(element % 2 == 0) swap_num = element/2;
            if(element == 0) continue;
            
            //routeをスワップ
            count = 0;
            for(p = point[0]; p < point[0] + swap_num; p++){
                tmp = route[p];
                route[p] = route[ point[1]-count ];
                route[ point[1]-count ] = tmp;
                count++;
            }
            
            //ルートのコストを計算
            result = Get_Circuit_Cost(distance, route, x, y);
            
            //改悪解に移動
            D = result - init;
            random = rand()%1001;
            random = (float)random/100000 + 0.99;
            if(D <= 0 || ( D > 0 && exp(-1*D/T) >= random ) ){
                flag = true;
                min = result;
                init = result;
                break;
            }
            //パラメータ更新
            T = ALFA * T;
            T = (int)T;
            if(T > INIT_T/2) {
                R = INIT_R * ( (1.5*INIT_T - T) / INIT_T ) * ( (1.5*INIT_T - T) / INIT_T ) * ( (1.5*INIT_T - T) / INIT_T ) * ( (1.5*INIT_T - T) / INIT_T );
            }
            if(T < INIT_T/2) {
                R = INIT_R * ( (0.5*INIT_T + T)  / INIT_T ) * ( (0.5*INIT_T + T)  / INIT_T ) * ( (0.5*INIT_T + T)  / INIT_T ) * ( (0.5*INIT_T + T)  / INIT_T );
            }
            if(T == INIT_T/2) R = INIT_R;
            R = (int)R;
            
            //route戻す
            count = 0;
            for(p = point[0]; p < point[0] + swap_num; p++){
                tmp = route[p];
                route[p] = route[ point[1]-count ];
                route[ point[1]-count ] = tmp;
                count++;
            }
        }
    }
    //計算時間計測終わり
    t2 = getrusage_sec();
    execution_time = t2 - t1;
   
    return min;
}


long Get_Distance(double x1, double y1, double x2, double y2){
    long dis_x = x1 - x2;
    long dis_y = y1 - y2;
    return (long)sqrt( dis_x*dis_x + dis_y*dis_y );
}


//ルートを取得
void Get_Route_Array(int route[n]){
    int i, seed;
    
    //ルート表の初期化
    for(i = 0; i < n; i++) route[i] = 0;
    
    //ランダムなルートを設定
    set_seed++;
    sgenrand(set_seed);
    for(i = 0; i < n; i++){
        double rand = genrand()*1000;
        rand = (int)rand % n;
        while(route[(int)rand] != 0) {
            rand = ((int)rand + 1)%n;
        }
        route[(int)rand] = i;
    }
}


//ルートのコストを計算
double Get_Circuit_Cost(long **distance, int route[], double x[], double y[]){
    double cost = 0;
    int i;
    
    for(i = 0; i < n; i++){
        cost += distance[ route[i] ][ route[(i+1)%n] ];
    }
    cost += distance[ route[0] ][ route[n-1] ];
    
    try_num++;
    if(MODE != 3) GNU_Plot(try_num, route, x, y);
    
    return cost;
}


//プロット
void GNU_Plot(int interation, int route[], double x[], double y[]){
    FILE *plot_fp;
    char filename[50];
    int i,j;
    
    
    if(interation % PLOT_INTERVAL == 0){
        if(MODE == 1) sprintf(filename, "random-%d.dat", interation);
        if(MODE == 2) sprintf(filename, "hill-climbing_%s-%d.dat", vicinity, interation);
        if(MODE == 3) sprintf(filename, "simulated-annealing-%d.dat", interation);
        printf("%s\n", filename);
        if( (plot_fp = fopen(filename, "w") ) == NULL){
            printf("File Open Error: %s\n", filename);
            exit(1);
        }
        fprintf(plot_fp, "#問題：%s\n", data_file);
        fprintf(plot_fp, "#解法：%s\n", solve_method);
        if(MODE == 2 || MODE == 3) fprintf(plot_fp, "#近傍：%s\n", vicinity);
        if(MODE == 3) {
            fprintf(plot_fp, "#初期温度：%d\n", INIT_T);
            fprintf(plot_fp, "#初期反復回数：%d\n", INIT_R);
            fprintf(plot_fp, "#温度減少率：%.3f\n", ALFA);
            fprintf(plot_fp, "#暫定解：%d\n", SA_ans);
        }
        if(MODE != 3) fprintf(plot_fp, "#繰り返し回数：%d\n", interation);
        if(MODE == 3) fprintf(plot_fp, "#実行時間：%.3f\n", execution_time);
        for(j = 0; j < 51; j++){
            fprintf(plot_fp, "%lf %lf\n", x[route[j]], y[route[j]]);
        }
        fclose(plot_fp);
    }
    
}


//計測時間測定
double getrusage_sec(){
    struct rusage t;
    struct timeval tv;
    getrusage(RUSAGE_SELF, &t);
    tv = t.ru_utime;
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}
