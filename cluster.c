// File: cluster.c
// Subject: IZP
// Project: #2
// Author: Andrii Klymenko, FIT VUT
// Login: xklyme00
// Date: 22.7.2023

/**
 * Kostra programu pro 2. projekt IZP 2022/23
 *
 * Jednoducha shlukova analyza: 2D nejblizsi soused.
 * Single linkage
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> // sqrtf
#include <stdbool.h>
#include <string.h>
#include <time.h>

/*****************************************************************
 * Ladici makra. Vypnout jejich efekt lze definici makra
 * NDEBUG, napr.:
 *   a) pri prekladu argumentem prekladaci -DNDEBUG
 *   b) v souboru (na radek pred #include <assert.h>
 *      #define NDEBUG
 */
#ifdef NDEBUG
#define debug(s)
#define dfmt(s, ...)
#define dint(i)
#define dfloat(f)
#else

// vypise ladici retezec
#define debug(s) printf("- %s\n", s)

// vypise formatovany ladici vystup - pouziti podobne jako printf
#define dfmt(s, ...) printf(" - "__FILE__":%u: "s"\n",__LINE__,__VA_ARGS__)

// vypise ladici informaci o promenne - pouziti dint(identifikator_promenne)
#define dint(i) printf(" - " __FILE__ ":%u: " #i " = %d\n", __LINE__, i)

// vypise ladici informaci o promenne typu float - pouziti
// dfloat(identifikator_promenne)
#define dfloat(f) printf(" - " __FILE__ ":%u: " #f " = %g\n", __LINE__, f)

#endif

#define MAX_CLUSTER_NUMBER 10000  // maximum number of clusters that can be processed by the program
#define MAX_LINE_LENGTH 15  // in the case when the line is "10000 1000 1000"
#define MAX_LINE_BUFFER_LENGTH MAX_LINE_LENGTH + 2  // + '\n' + '\0'
#define DELIMITER_CHAR ' '  // space is a delimiter
#define DELIMITER_STRING " "  // space is a delimiter (for strtok function)
#define MAX_CLUSTER_DISTANCE 1414.21356 + 1.0  // in the case when first objects has [0, 0] coordinates and second objects has [1000, 1000] coordinates + 1
#define MIN_CLUSTER_DISTANCE 0.0 - 1.0  // in the case when two objects have the same coordinates - 1

/*****************************************************************
 * Deklarace potrebnych datovych typu:
 *
 * TYTO DEKLARACE NEMENTE
 *
 *   struct obj_t - struktura objektu: identifikator a souradnice
 *   struct cluster_t - shluk objektu:
 *      pocet objektu ve shluku,
 *      kapacita shluku (pocet objektu, pro ktere je rezervovano
 *          misto v poli),
 *      ukazatel na pole shluku.
 */

// structure that contains an information about the object
typedef struct obj_t {
    int id; // object id
    float x; // object x-coordinate
    float y; // object y-coordinate
} obj_t;

typedef struct cluster_t {
    int size; // number of the objects in the cluster
    int capacity; // number of the objects that cluster can have without additional memory allocation
    struct obj_t *obj; // an array of the cluster objects
} cluster_t;

// structure that contains program arguments
typedef struct arguments_t {
    char *filename;  // file which contains objects
    char flag;  // specifies which clustering algorithm will be used
    int required_clusters;  // final number of clusters
} arguments_t;

// a pointer to the function that calculates a distance between two objects depending on the clustering algorithm
// it has two parameters (cluster_t * and cluster_t *) and returns a float value
typedef float (*distanceFunction)(cluster_t *, cluster_t *);

/*****************************************************************
 * Deklarace potrebnych funkci.
 *
 * PROTOTYPY FUNKCI NEMENTE
 *
 * IMPLEMENTUJTE POUZE FUNKCE NA MISTECH OZNACENYCH 'TODO'
 *
 */

/*
 Inicializace shluku 'c'. Alokuje pamet pro cap objektu (kapacitu).
 Ukazatel NULL u pole objektu znamena kapacitu 0.
*/

void *init_cluster(cluster_t *c, int cap)
{
    assert(c != NULL);
    assert(cap >= 0);

    if(cap > 0)
    {
        c->obj = (obj_t *) malloc(sizeof(obj_t) * cap);

        if(c->obj == NULL)
            return NULL;
    }
    else
        c->obj = NULL;

    c->capacity = cap;
    c->size = 0;

    return c;
}

/*
 Odstraneni vsech objektu shluku a inicializace na prazdny shluk.
 */
void clear_cluster(cluster_t *c)
{
    if(c->obj != NULL)
        free(c->obj);

    init_cluster(c, 0);
}

/// Chunk of cluster objects. Value recommended for reallocation.
const int CLUSTER_CHUNK = 10;

/*
 Zmena kapacity shluku 'c' na kapacitu 'new_cap'.
 */
cluster_t *resize_cluster(cluster_t *c, int new_cap)
{
    // TUTO FUNKCI NEMENTE
    assert(c);
    assert(c->capacity >= 0);
    assert(new_cap >= 0);

    if (c->capacity >= new_cap)
        return c;

    size_t size = sizeof(obj_t) * new_cap;

    void *arr = realloc(c->obj, size);
    if (arr == NULL)
        return NULL;

    c->obj = (obj_t*)arr;
    c->capacity = new_cap;
    return c;
}

/*
 Prida objekt 'obj' na konec shluku 'c'. Rozsiri shluk, pokud se do nej objekt
 nevejde.
 */

void *append_cluster(cluster_t *c, obj_t obj)
{
    if(c->capacity == c->size)
    {
        if(resize_cluster(c, c->capacity + CLUSTER_CHUNK) == NULL)
            return NULL;
    }

    c->obj[c->size++] = obj;
    return c;
}

/*
 Seradi objekty ve shluku 'c' vzestupne podle jejich identifikacniho cisla.
 */
void sort_cluster(cluster_t *c);

/*
 Do shluku 'c1' prida objekty 'c2'. Shluk 'c1' bude v pripade nutnosti rozsiren.
 Objekty ve shluku 'c1' budou serazeny vzestupne podle identifikacniho cisla.
 Shluk 'c2' bude nezmenen.
 */

void *merge_clusters(cluster_t *c1, cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c2 != NULL);

    for(int j = 0; j < c2->size; j++)
        append_cluster(c1, c2->obj[j]);

    sort_cluster(c1);
    return c1;
}

/**********************************************************************/
/* Prace s polem shluku */

/*
 Odstrani shluk z pole shluku 'carr'. Pole shluku obsahuje 'narr' polozek
 (shluku). Shluk pro odstraneni se nachazi na indexu 'idx'. Funkce vraci novy
 pocet shluku v poli.
*/
int remove_cluster(cluster_t *carr, int narr, int idx)
{
    assert(idx < narr);
    assert(narr > 0);

    clear_cluster(&carr[idx]);
    narr--;

    for(int i = idx; i < narr; i++)
        carr[i] = carr[i + 1];

    init_cluster(&carr[narr], 0);

    return narr;
}

/*
 Pocita Euklidovskou vzdalenost mezi dvema objekty.
 */
float obj_distance(obj_t *o1, obj_t *o2)
{
    assert(o1 != NULL);
    assert(o2 != NULL);

    return sqrt((o1->x - o2->x) * (o1->x - o2->x) + (o1->y - o2->y) * (o1->y - o2->y));
}

/*
 Pocita vzdalenost dvou shluku.
*/
// single linkage
// the distance of two clusters is equal to the smallest distance of any two objects from both clusters
float cluster_distance_single(cluster_t *c1, cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float min = MAX_CLUSTER_DISTANCE;

    for(int i = 0; i < c1->size; i++)
    {
        for(int j = 0; j < c2->size; j++)
        {
            float distance = obj_distance(&c1->obj[i], &c2->obj[j]);

            if(distance < min)
                min = distance;
        }
    }

    return min;
}

// complete linkage
// the distance of two clusters is equal to the largest distance of any two objects from both clusters
float cluster_distance_complete(cluster_t *c1, cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float max = MIN_CLUSTER_DISTANCE;

    for(int i = 0; i < c1->size; i++)
    {
        for(int j = 0; j < c2->size; j++)
        {
            float distance = obj_distance(&c1->obj[i], &c2->obj[j]);

            if(distance > max)
                max = distance;
        }
    }

    return max;
}

// average linkage
// the distance between each pair of objects in each cluster are added up and divided by the number of pairs to get an average inter-cluster distance
float cluster_distance_average(cluster_t *c1, cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    float avg = MIN_CLUSTER_DISTANCE;

    for(int i = 0; i < c1->size; i++)
        for(int j = 0; j < c2->size; j++)
            avg += obj_distance(&c1->obj[i], &c2->obj[j]);

    return avg / (float) (c1->size * c2->size);
}

/*
 Funkce najde dva nejblizsi shluky. V poli shluku 'carr' o velikosti 'narr'
 hleda dva nejblizsi shluky. Nalezene shluky identifikuje jejich indexy v poli
 'carr'. Funkce nalezene shluky (indexy do pole 'carr') uklada do pameti na
 adresu 'c1' resp. 'c2'.
*/
void find_neighbours(cluster_t *carr, int narr, int *c1, int *c2, distanceFunction get_distance)
{
    assert(narr > 0);

    float min = MAX_CLUSTER_DISTANCE;

    for(int i = 0; i < narr - 1; i++)
    {
        for(int j = i + 1; j < narr; j++)
        {
            float distance = get_distance(&carr[i], &carr[j]);

            if(distance < min)
            {
                min = distance;
                *c1 = i;
                *c2 = j;
            }
        }
    }
}

// pomocna funkce pro razeni shluku
static int obj_sort_compar(const void *a, const void *b)
{
    // TUTO FUNKCI NEMENTE
    const obj_t *o1 = (const obj_t *)a;
    const obj_t *o2 = (const obj_t *)b;
    if (o1->id < o2->id) return -1;
    if (o1->id > o2->id) return 1;
    return 0;
}

/*
 Razeni objektu ve shluku vzestupne podle jejich identifikatoru.
*/
void sort_cluster(cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    qsort(c->obj, c->size, sizeof(obj_t), &obj_sort_compar);
}

/*
 Tisk shluku 'c' na stdout.
*/
void print_cluster(cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    for (int i = 0; i < c->size; i++)
    {
        if (i) putchar(' ');
        printf("%d[%g,%g]", c->obj[i].id, c->obj[i].x, c->obj[i].y);
    }
    putchar('\n');
}

// checks if string contain an integer from the interval [1, 10000] inclusively
bool checkNumber(char *str, int *number)
{
    char *end_ptr;
    long result = strtol(str, &end_ptr, 10);

    if(*end_ptr != '\0' || result < 1 || result > MAX_CLUSTER_NUMBER)
        return false;

    *number = (int) result;
    return true;
}

// checks if the line contains a '\n' character
// i.e. it contains at most 15 characters ("10000 1000 1000")
bool checkLineLength(char *line)
{
    if(strchr(line, '\n') == NULL)
        return false;

    line[strcspn(line, "\n")] = '\0'; // replaces '\n' with '\0'
    return true;
}

// checks first line in the file
bool checkFirstLine(char *line, int *object_number)
{
    if(!checkLineLength(line))
    {
        fprintf(stderr, "Error! Invalid line length of the first line in the file\n");
        return false;
    }

    if(strncmp(line, "count=", 6) != 0)
    {
        fprintf(stderr, "Error! First line must start from the 'count='\n");
        return false;
    }

    char *number = line + 6;

    if(checkNumber(number, object_number))
        return true;

    fprintf(stderr, "Error! Number of the objects in the file must be a positive integer\n");
    return false;
}

// checks if there are two delimiters (spaces) following each other
bool checkDelimiterSequence(char *line)
{
    for(int i = 0; line[i] != '\0'; i++)
        if(line[i] == DELIMITER_CHAR && line[i + 1] == DELIMITER_CHAR)
            return false;

    return true;
}

// check if line contains only two delimiters (spaces)
// id x y
bool checkDelimiterNumber(char *line)
{
    int delimiter_number = 0;

    for(int i = 0; line[i] != '\0'; i++)
    {
        if(line[i] == DELIMITER_CHAR)
        {
            delimiter_number++;

            if(delimiter_number > 2)
                return false;
        }
    }

    return true;
}

// checks if the id is unique in the cluster/object array
bool isIdUnique(cluster_t *cluster_arr, obj_t *object_arr, int size, int id, char flag)
{
    if(flag != 'k')
    {
        for(int i = 0; i < size - 1; i++)  // we know that each cluster has exactly one object
            if(cluster_arr[i].obj->id == id)
                return false;

        return true;
    }

    // else

    for(int i = 0; i < size - 1; i++)
        if(object_arr[i].id == id)
            return false;

    return true;
}

// checks if string contains a float value with a zero decimal which is from the interval [0, 1000] inclusively
bool checkCoordinate(char *str, float *num)
{
    char *end_ptr;
    *num = strtof(str, &end_ptr);

    double fractional_part = modf(*num, &fractional_part);

    return fractional_part == 0.0 && *end_ptr == '\0' && *num >= 0.0 && *num <= 1000.0;
}

// checks line declaring an object in the file
bool checkObjectLine(char *line, cluster_t *cluster_arr, obj_t *object_arr, int line_cnt, char flag)
{
    if(!checkLineLength(line))
    {
        fprintf(stderr, "Error! Invalid line length in the line no. %d declaring an object'\n", line_cnt + 1);
        return false;
    }

    if(!checkDelimiterNumber(line) || !checkDelimiterSequence(line))
    {
        fprintf(stderr, "Error! There should be exactly two not following each other delimiters in the line declaring an object\n");
        return false;
    }

    obj_t *obj;

    // if program performs k-means clustering the object from the file will be copied to the array of the objects
    // otherwise it will be copied to the array of clusters (and it will be the only cluster object for a while)
    if(flag != 'k')
        obj = cluster_arr[line_cnt - 1].obj;
    else
        obj = &(object_arr[line_cnt - 1]);

    char *token = strtok(line, DELIMITER_STRING); // object id (string)

    if(!checkNumber(token, &(obj->id)))
    {
        fprintf(stderr, "Error! Invalid object id on the line no. %d", line_cnt + 1);
        return false;
    }

    if(!isIdUnique(cluster_arr, object_arr, line_cnt, obj->id, flag))
    {
        fprintf(stderr, "Error! Every object id must be unique\n");
        return false;
    }

    token = strtok(NULL, DELIMITER_STRING); // object x-coordinate (string)

    if(!checkCoordinate(token, &(obj->x)))
    {
        fprintf(stderr, "Error! Invalid object x coordinate on the line no. %d", line_cnt + 1);
        return false;
    }

    token = strtok(NULL, DELIMITER_STRING); // object y-coordinate (string)

    if(checkCoordinate(token, &(obj->y)))
    {
        if(flag != 'k')
            (cluster_arr[line_cnt - 1].size)++; // increment the number of the objects in the cluster

        return true;
    }

    fprintf(stderr, "Error! Invalid object y coordinate on the line no. %d", line_cnt + 1);
    return false;
}

// frees a memory that was dynamically allocated
void destroy(cluster_t *cluster_arr, int cluster_arr_size, obj_t *object_arr)
{
    if(cluster_arr != NULL)
    {
        for(int i = 0; i < cluster_arr_size; i++)
            clear_cluster(&cluster_arr[i]);

        free(cluster_arr);
        cluster_arr = NULL;
    }

    if(object_arr != NULL)
    {
        free(object_arr);
        object_arr = NULL;
    }
}

/*
 Tisk pole shluku. Parametr 'carr' je ukazatel na prvni polozku (shluk).
 Tiskne se prvnich 'narr' shluku.
*/
void print_clusters(cluster_t *carr, int narr)
{
    printf("Clusters:\n");
    for (int i = 0; i < narr; i++)
    {
        printf("cluster %d: ", i);
        print_cluster(&carr[i]);
    }
}

bool initAllClusters(cluster_t *cluster_arr, int cluster_arr_size)
{
    // firstly initialize a cluster to {.size = 0, .cap = 0, .obj = NULL}
    // in order to prevent freeing memory that wasn't allocated
    for(int i = 0; i < cluster_arr_size; i++)
        init_cluster(&cluster_arr[i], 0);

    // then initialize every cluster with capacity 1
    for(int i = 0; i < cluster_arr_size; i++)
        if(init_cluster(&cluster_arr[i], 1) == NULL)
            return false;

    return true;
}

bool init(cluster_t **cluster_arr, obj_t **object_arr, int arr_size, arguments_t *a)
{
    if(a->flag != 'k')
    {
        // allocate memory for an array of the clusters and initialize all clusters
        *cluster_arr = (cluster_t *) malloc(arr_size * sizeof(cluster_t));
        return *cluster_arr != NULL && initAllClusters(*cluster_arr, arr_size);
    }

    // if flag == 'k'
    // allocate memory for an array of the objects and clusters and initialize all clusters
    *object_arr = (obj_t *) malloc(arr_size * sizeof(obj_t));
    *cluster_arr = (cluster_t *) malloc(a->required_clusters * sizeof(cluster_t));

    return *object_arr != NULL && *cluster_arr != NULL && initAllClusters(*cluster_arr, a->required_clusters);
}

bool processFile(cluster_t **cluster_arr, obj_t **object_arr, int *arr_size, FILE *f, arguments_t *a)
{
    char line[MAX_LINE_BUFFER_LENGTH]; // string that contains a line from the file

    if((fgets(line, MAX_LINE_BUFFER_LENGTH, f)) != NULL)  // read first line
    {
        if(!checkFirstLine(line, arr_size))
            return false;

        if(!init(cluster_arr, object_arr, *arr_size, a))
        {
            if(a->flag != 'k')
                fprintf(stderr, "Error! Couldn't allocate memory for an array of clusters/initialize a cluster\n");
            else
                fprintf(stderr, "Error! Couldn't allocate memory for an array of objects/clusters/initialize a cluster\n");

            return false;
        }
    }
    else
    { // if fgets function returned NULL
        fprintf(stderr, "Error! File '%s' is empty\n", a->filename);
        return false;
    }

    int line_cnt = 1; // number of the lines in the file

    while((fgets(line, MAX_LINE_BUFFER_LENGTH, f)) != NULL)
    {
        // if the program already read specified number of the objects stop reading file
        if(line_cnt == *arr_size + 1)  // line_cnt == number of objects + 1 (first line 'count=x')
            break;

        if(!checkObjectLine(line, *cluster_arr, *object_arr, line_cnt, a->flag))
            return false;

        line_cnt++;
    }

    // if there are fewer objects than first line specified
    if(line_cnt != *arr_size + 1)
    {
        fprintf(stderr, "Error! Expected %d objects in the file '%s', but got %d\n", *arr_size, a->filename, line_cnt - 1);
        return false;
    }

    fclose(f);
    return true;
}

/*
 Ze souboru 'filename' nacte objekty. Pro kazdy objekt vytvori shluk a ulozi
 jej do pole shluku. Alokuje prostor pro pole vsech shluku a ukazatel na prvni
 polozku pole (ukalazatel na prvni shluk v alokovanem poli) ulozi do pameti,
 kam se odkazuje parametr 'arr'. Funkce vraci pocet nactenych objektu (shluku).
 V pripade nejake chyby uklada do pameti, kam se odkazuje 'arr', hodnotu NULL.
*/

int load_clusters(cluster_t **cluster_arr, obj_t **object_arr, arguments_t *a)
{
    assert(cluster_arr != NULL);

    FILE *f = fopen(a->filename, "r");

    if(f == NULL)
    {
        fprintf(stderr, "Error! Couldn't open a file '%s'\n", a->filename);
        return -1;
    }

    // number of objects in the object array/number of clusters in cluster array (depends on clustering algorithm)
    int arr_size = 0;

    if(!processFile(cluster_arr, object_arr, &arr_size, f, a))
    {
        // if some error occurred and the program was performing k-means algorithm, only 'required_clusters' clusters
        // need to be freed
        if(a->flag == 'k')
            arr_size = a->required_clusters;

        destroy(*cluster_arr, arr_size, *object_arr);
        fclose(f);
        return -1;
    }

    return arr_size;
}

// check if string contains a valid program flag
// list of the valid flags:
// '-c' - complete linkage
// '-a' - average linkage
// '-k' - k-means
// '-s' - single linkage
bool checkFlag(char *str, char *flag)
{
    if(strcmp(str, "-c") == 0 || strcmp(str, "-a") == 0 || strcmp(str, "-k") == 0 || strcmp(str, "-s") == 0)
    {
        *flag = str[1];
        return true;
    }

    return false;
}

// ./executable filename [N] [flag]
// parses program arguments
// returns false if invalid argument was encountered
// returns true if all arguments were valid
bool parseArguments(int argc, char *argv[], arguments_t *a)
{
    if(argc < 2 || argc > 4)
    {
        fprintf(stderr, "Error! Invalid number of the program arguments\n");
        return false;
    }

    a->filename = argv[1];

    //  ./executable filename
    if(argc == 2)
        return true;

    // ./executable filename N or ./executable filename flag
    if(argc == 3)
    {
        if(checkNumber(argv[2], &a->required_clusters))
            return true;

        if(checkFlag(argv[2], &a->flag))
            return true;

        fprintf(stderr, "Error! Invalid third program argument '%s'\n", argv[2]);
        return false;
    }

    // if argc == 4
    // ./executable filename N flag
    if(!checkNumber(argv[2], &a->required_clusters))
    {
        fprintf(stderr, "Error! Invalid third program argument '%s'\n", argv[2]);
        return false;
    }

    if(checkFlag(argv[3], &a->flag))
        return true;

    fprintf(stderr, "Error! Invalid fourth program argument '%s'\n", argv[3]);
    return false;
}

// checks if an array of integers contains a specific value
bool containDuplicate(int *arr, int size, int value)
{
    for(int i = 0; i < size; i++)
        if(arr[i] == value)
            return true;

    return false;
}

// on a success, returns an array with random numbers from the range [1, objects_arr_size - 1] inclusively
// otherwise returns NULL
int *getRandomNumbers(int arr_size, int object_arr_size)
{
    int *random_numbers = (int *) malloc(sizeof(int) * arr_size); // allocate memory

    if(random_numbers == NULL) // check if allocation was successful
    {
        fprintf(stderr, "Error! Couldn't allocate memory for an array of centroids ids\n");
        return NULL;
    }

    srand(time(NULL)); // set the seed to get truly random numbers each time the program is run

    for(int i = 0; i < arr_size; )
    {
        // get a random number from the range [1, object_arr_size - 1] inclusively (getting and index of object array)
        int random_number = rand() % object_arr_size;

        // check if the generated random number wasn't generated before, i.e. it is unique
        if(!containDuplicate(random_numbers, i, random_number))
        {
            random_numbers[i] = random_number; // write a value to random_numbers arr
            i++;
        }
    }

    return random_numbers; // return an array with random numbers
}

// initializes centroid array
bool initializeCentroids(obj_t **centroid_arr, int centroid_arr_size, obj_t *object_arr, int object_arr_size)
{
    *centroid_arr = (obj_t *) malloc(sizeof(obj_t) * centroid_arr_size); // allocate memory for centroid arr

    if(*centroid_arr == NULL) // check if allocation was successfully
    {
        fprintf(stderr, "Error! Couldn't allocate memory for an array of centroids\n");
        return false;
    }

    // get an array with indexes of the objects from the object array that will become centroids
    int *indexes = getRandomNumbers(centroid_arr_size, object_arr_size);

    if(indexes == NULL) // couldn't allocate a memory for an array of random numbers
    {
        fprintf(stderr, "Error! Couldn't allocate memory for an array of random numbers\n");
        free(centroid_arr); // free memory that was allocated for a centroid array
        return false;
    }

    // copy objects that have to become centroids to centroid array
    for(int i = 0; i < centroid_arr_size; i++)
        (*centroid_arr)[i] = object_arr[indexes[i]];

    free(indexes); // free memory that was allocated for an array of random numbers
    return true;
}

// calculates the distance between an object and every centroid from centroid array
// returns the smallest distance
float getMinDistance(obj_t *obj, obj_t *centroid_arr, int centroid_arr_size)
{
    float min = MAX_CLUSTER_DISTANCE;

    for(int i = 0; i < centroid_arr_size; i++)
    {
        float distance = obj_distance(obj, &centroid_arr[i]);

        if(distance < min)
            min = distance;
    }

    return min;
}

// assigns an object to the cluster based on the distance to the cluster centroid
bool assignObjectToCluster(obj_t *obj, float distance, obj_t *centroid_arr, int centroid_arr_size, cluster_t *cluster_arr)
{
    for(int i = 0; i < centroid_arr_size; i++)
    {
        if(distance == obj_distance(obj, &centroid_arr[i]))
        {
            if(append_cluster(&cluster_arr[i], *obj) != NULL)
                return true;

            break;
        }
    }

    return false;
}

// shifts all objects that are to the right of this index 'idx' and decrements cluster size
void shiftCluster(cluster_t *cluster, int idx)
{
    for(int i = idx; i < cluster->size - 1; i++)
        cluster->obj[i] = cluster->obj[i + 1];

    (cluster->size)--;
}

// removes an object from the array of clusters based on object id
void removeObjectFromClusterArr(cluster_t *cluster_arr, int cluster_arr_size, obj_t *obj)
{
    for(int i = 0; i < cluster_arr_size; i++)
    {
        for(int j = 0; j < cluster_arr[i].size; j++)
            if(obj->id == cluster_arr[i].obj[j].id)
                shiftCluster(&cluster_arr[i], j);
    }
}

// assigns every object to the cluster
bool assignObjectsToClusters(obj_t *object_arr, int object_arr_size, cluster_t *cluster_arr, int required_clusters, obj_t *centroid_arr)
{
    for(int i = 0; i < object_arr_size; i++)
    {
        float min = getMinDistance(&object_arr[i], centroid_arr, required_clusters);

        obj_t obj = object_arr[i];

        // remove the object from the previous cluster
        removeObjectFromClusterArr(cluster_arr, required_clusters, &obj);

        // assign the object to the new cluster
        if(!assignObjectToCluster(&obj, min, centroid_arr, required_clusters, cluster_arr))
        {
            free(centroid_arr);
            return false;
        }
    }

    return true;
}

// updates cluster centroid
// if centroid wasn't updates returns false, otherwise returns true
bool updateClusterCentroid(obj_t *centroid, cluster_t *cluster)
{
    obj_t tmp = *centroid;
    centroid->x = centroid->y = 0.0;

    for(int i = 0; i < cluster->size; i++)
    {
        centroid->x += cluster->obj[i].x;
        centroid->y += cluster->obj[i].y;
    }

    centroid->x /= (float) cluster->size;
    centroid->y /= (float) cluster->size;

    if(tmp.x != centroid->x || tmp.y != centroid->y)
        return true;

    return false;
}

// updates all cluster centroids
bool updateClusterCentroids(obj_t *centroid_arr, int centroid_arr_size, cluster_t *cluster_arr)
{
    for(int i = 0; i < centroid_arr_size; i++)
        if(!updateClusterCentroid(&centroid_arr[i], &cluster_arr[i]))
            return false;

    return true;
}

// implementation of k-means clustering algorithm
bool kMeansClustering(obj_t *object_arr, int object_arr_size, cluster_t *cluster_arr, int required_clusters)
{
    obj_t *centroid_arr = NULL; // an array of centroids

    // initialize centroids and make a first assignment of the objects to the clusters
    if(!initializeCentroids(&centroid_arr, required_clusters, object_arr, object_arr_size) ||
       !assignObjectsToClusters(object_arr, object_arr_size, cluster_arr, required_clusters, centroid_arr))
        return false;

    // reassign objects to centroids until centroids don't change
    while(updateClusterCentroids(centroid_arr, required_clusters, cluster_arr))
    {
        if(!assignObjectsToClusters(object_arr, object_arr_size, cluster_arr, required_clusters, centroid_arr))
            return false;
    }

    free(centroid_arr);
    return true;
}

// implementation of single/complete/average linkage clustering algorithms
bool defaultClustering(int *arr_size, int required_clusters, cluster_t *cluster_arr, distanceFunction get_distance)
{
    while(*arr_size != required_clusters)
    {
        int c1, c2; // indexes of the clusters in the cluster array

        // get the indexes of the clusters that have to be merged
        find_neighbours(cluster_arr, *arr_size, &c1, &c2, get_distance);

        // cluster with index c2 will be merged to the cluster with index c1
        if(merge_clusters(&cluster_arr[c1], &cluster_arr[c2]) == NULL)
        {
            fprintf(stderr, "Error! Couldn't merge two clusters\n");
            return false;
        }

        // after merging remove cluster with index c2
        *arr_size = remove_cluster(cluster_arr, *arr_size, c2);
    }

    return true;
}

// gets required number of clusters
int finalClustering(int *arr_size, cluster_t *cluster_arr, distanceFunction get_distance, obj_t *object_arr, arguments_t *a)
{
    if(a->required_clusters > *arr_size)
    {
        fprintf(stderr, "Error! Third program argument %d is greater than number of the objects (%d)\n", a->required_clusters, *arr_size);
        return -1;
    }

    if(a->flag != 'k')
    {
        if(!defaultClustering(arr_size, a->required_clusters, cluster_arr, get_distance))
            return -1;
    }
    else
    {
        if(!kMeansClustering(object_arr, *arr_size, cluster_arr, a->required_clusters))
        {
            *arr_size = a->required_clusters; // need to free only 'required_clusters' clusters
            return -1;
        }

        *arr_size = a->required_clusters;
    }

    print_clusters(cluster_arr, *arr_size);
    return 0;
}

int main(int argc, char *argv[])
{
    cluster_t *cluster_arr = NULL; // an array of clusters
    obj_t *object_arr = NULL; // an array of objects (for k-means clustering)

    arguments_t a = {.required_clusters = 1, .flag = 's'}; // program arguments

    if(!parseArguments(argc, argv, &a))
        return -1;

    // in the case when program performs k-means clustering, arr_size represents number of the objects
    // in the object_arr
    // otherwise it represents number of the clusters
    int arr_size = load_clusters( &cluster_arr, &object_arr, &a);

    if(arr_size == -1)
        return -1;

    distanceFunction get_distance = NULL;

    if(a.flag == 's') // single linkage
        get_distance = cluster_distance_single;
    else if(a.flag == 'c') // complete linkage
        get_distance = cluster_distance_complete;
    else if(a.flag == 'a') // average linkage
        get_distance = cluster_distance_average;

    int result = finalClustering(&arr_size, cluster_arr, get_distance, object_arr, &a);

    destroy(cluster_arr, arr_size, object_arr);
    return result;
}
