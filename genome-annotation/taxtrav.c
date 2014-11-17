/*
Copyright (c) 2013-2014, Daniel S. Standage <daniel.standage@gmail.com>

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

// For compiling try: gcc -Wall -O3 -o taxtrav taxtrav.c

#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//------------------------------------------------------------------------------
// Basic hash map data structure
//------------------------------------------------------------------------------
typedef struct _hashmap_item_
{
  char *key;
  void *value;
  struct _hashmap_item_ *next;
} HashmapItem;
typedef struct
{
  int size;
  int usage;
  HashmapItem **table;
} Hashmap;

// Method prototypes
int hashmap_add(Hashmap *hash, char *key, void *value);
void hashmap_delete(Hashmap *hash, void (*valuefreefunc)(void *));
void *hashmap_get(Hashmap *hash, char *key);
unsigned int hashmap_hash(Hashmap *hash, char *key);
char **hashmap_keys(Hashmap *hash);
Hashmap *hashmap_new(int size);
int hashmap_size(Hashmap *hash);

// Method implementations
int hashmap_add(Hashmap *hash, char *key, void *value)
{
  HashmapItem *item;
  unsigned int hashval;

  if(hashmap_get(hash, key) != NULL)
    return 0;

  hashval = hashmap_hash(hash, key);
  item = malloc(sizeof(HashmapItem));
  if(item == NULL)
    return 0;
  item->key = strdup(key);
  item->value = value;
  item->next = hash->table[hashval];
  hash->table[hashval] = item;
  hash->usage += 1;

  return 1;
}
void hashmap_delete(Hashmap *hash, void (*valuefreefunc)(void *))
{
  HashmapItem *item, *temp;
  int i;

  if(hash == NULL)
    return;

  for(i = 0; i < hash->size; i++)
  {
    item = hash->table[i];
    while(item != NULL)
    {
      temp = item;
      item = item->next;
      free(temp->key);
      if(valuefreefunc != NULL)
        valuefreefunc(temp->value);
      free(temp);
    }
  }

  free(hash->table);
  free(hash);
}
void *hashmap_get(Hashmap *hash, char *key)
{
  HashmapItem *item;
  unsigned int hashval = hashmap_hash(hash, key);
  for(item = hash->table[hashval]; item != NULL; item = item->next)
  {
    if(strcmp(key, item->key) == 0)
      return item->value;
  }
  return NULL;
}
unsigned int hashmap_hash(Hashmap *hash, char *key)
{
  unsigned int hashval = 0;
  for(; *key != '\0'; key++)
    hashval = *key + (hashval << 5) - hashval;
  return hashval % hash->size;
}
char **hashmap_keys(Hashmap *hash)
{
  char **keys;
  HashmapItem *item;
  int i, j;

  keys = malloc(sizeof(char *) * hash->usage);
  j = 0;
  for(i = 0; i < hash->size; i++)
  {
    item = hash->table[i];
    while(item != NULL)
    {
      keys[j++] = item->key;
      item = item->next;
    }
  }

  return keys;
}
Hashmap *hashmap_new(int size)
{
  Hashmap *hash;
  int i;

  if(size < 1)
    return NULL;

  hash = malloc(sizeof(HashmapItem));
  if(hash == NULL)
    return NULL;

  hash->table = malloc(sizeof(HashmapItem *) * size);
  if(hash->table == NULL)
    return NULL;

  for(i = 0; i < size; i++)
    hash->table[i] = NULL;

  hash->size = size;
  hash->usage = 0;
  return hash;
}
int hashmap_size(Hashmap *hash)
{
  return hash->usage;
}


//------------------------------------------------------------------------------
// Data structure to store info about particular taxonomic units
//------------------------------------------------------------------------------
typedef struct _taxon_
{
  char *id;
  char *rank;
  char *name;
  char *parentid;
  struct _taxon_ *parent;
} Taxon;

// Method prototypes
void taxon_delete(Taxon *t);
Taxon *taxon_new(const char *taxid, const char *parentid, const char *rank);

// Method implementations
void taxon_delete(Taxon *t)
{
  free(t->id);
  free(t->rank);
  free(t->parentid);
  if(t->name != NULL)
    free(t->name);  
  free(t);
  t = NULL;
}
Taxon *taxon_new(const char *taxid, const char *parentid, const char *rank)
{
  Taxon *t = malloc(sizeof(Taxon));
  t->id = strdup(taxid);
  t->rank = strdup(rank);
  t->parentid = strdup(parentid);
  t->name = NULL;
  t->parent = NULL;
  return t;
}


//------------------------------------------------------------------------------
// Main implementation
//------------------------------------------------------------------------------
void print_usage(FILE *outstream)
{
  fprintf( outstream, "\n"
"taxtrav: traverse the taxonomy to determine the classification of a given\n"
"         taxon at a particular taxonomic rank\n\n"
"Usage: taxtrav [options] input.file\n"
"  Options:\n"
"    -c|--complement:     if the --filter option is enabled, report only\n"
"                         results that do NOT match the specified taxon;\n"
"                         for example, if '-c -f order=Hymenoptera' is\n"
"                         specified, program will only report results for\n"
"                         input taxa that do not fall within the order\n"
"                         Hymenoptera\n"
"    -d|--delim: STRING   list of delimiters separating values in the input\n"
"                         file; default is ' ,\\t\\n' character\n"
"    -f|--filter: STRING  restrict results with a filter of the form\n"
"                         'rank=taxon_name'; for example, 'order=Hymenoptera'\n"
"                         will only report results for input taxon IDs that\n"
"                         fall within the order Hymenoptera; default is no\n"
"                         filter\n"
"    -h|--help            print this help message and exit\n"
"    -m|--names: FILE     path to file containing taxon names; default is\n"
"                         './names.dmp'; can be downloaded from\n"
"                         ftp://ftp.ncbi.nih.gov/pub/taxonomy/\n"
"    -n|--nodes: FILE     path to file containing taxonomy graph; default is\n"
"                         './nodes.dmp'; can be downloaded from\n"
"                         ftp://ftp.ncbi.nih.gov/pub/taxonomy/\n"
"    -r|--rank: STRING    desired taxonomic rank to report for each input\n"
"                         taxon; default is 'order'; to see all possible\n"
"                         taxonomic ranks, try the following command\n\n"
"                            cut -f 3 -d '|' < nodes.dmp \\\n"
"                              | perl -ne 's/^\\s*(.+)\\s*$/$1/;print $_,\"\\n\"' \\\n"
"                              | sort | uniq -c | sort -rn \n\n"
"The input file is expected to be a tab-, space-, or comma-delimited text\n"
"file, where the first value of each line is a GenBank taxonomy ID.\n\n"
"The output is 3 comma-delimited values for each line of input: the GenBank\n"
"ID of the input taxon, the GenBank ID of the taxon's classification at the\n"
"requested taxonomic rank, and the human-readable classification of the taxon\n"
"at the requested taxonomic rank.\n\n" );
}

FILE *myfopen(const char *filename, const char *mode)
{
  FILE *fh = fopen(filename, mode);
  if(fh == NULL)
  {
    fprintf(stderr, "error opening file '%s'\n", filename);
    exit(1);
  }
  return fh;
}

char *trim_str(char *str)
{
  while(isspace(str[0]))
    str += 1;
  int i = strlen(str) - 1;
  while(isspace(str[i]) && i >= 0)
  {
    str[i] = '\0';
    i--;
  }
  return str;
}

Taxon *get_parent_by_rank(Taxon *t, const char *rank)
{
  if(strcmp(t->id, "1") == 0) // Base case: found root of taxonomy
    return t;
  if(strcmp(t->rank, rank) == 0) // Found desired level of taxonomy
    return t;
  return get_parent_by_rank(t->parent, rank); // Recurse
}

int main(int argc, char **argv)
{
  int complement          = 0;
  const char *delims      = ", \t\n";
  const char *filterrank  = NULL;
  const char *filtertaxon = NULL;
  const char *nodesfile   = "nodes.dmp";
  const char *namesfile   = "names.dmp";
  const char *targetrank  = "order";

  // Parse options
  int opt = 0;
  int optindex = 0;
  const char *optstr = "cd:f:hn:m:r:";
  char *temparg;
  const struct option init_options[] =
  {
    { "complement", no_argument,       NULL, 'c' },
    { "delim",      required_argument, NULL, 'd' },
    { "filter",     required_argument, NULL, 'f' },
    { "help",       no_argument,       NULL, 'h' },
    { "names",      required_argument, NULL, 'm' },
    { "nodes",      required_argument, NULL, 'n' },
    { "rank",       required_argument, NULL, 'r' },
    { NULL,         no_argument,       NULL, 0 },
  };

  for( opt = getopt_long(argc, argv, optstr, init_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv, optstr, init_options, &optindex) )
  {
    switch(opt)
    {
      case 'c':
        complement = 1;
        break;

      case 'd':
        delims = optarg;
        break;

      case 'f':
        temparg = strdup(optarg);
        filterrank = strtok(optarg, "=");
        filtertaxon = strtok(NULL, "=");
        if(filterrank == NULL || filtertaxon == NULL)
        {
          fprintf(stderr, "error parsing filter '%s'\n", temparg);
          return 1;
        }
        free(temparg);
        break;

      case 'h':
        print_usage(stdout);
        return 0;
        break;

      case 'm':
        namesfile = optarg;
        break;

      case 'n':
        nodesfile = optarg;
        break;

      case 'r':
        targetrank = optarg;
        break;

      default:
        break;
    }
  }

  // Validate options and arguments
  int filenum = argc - optind;
  if(filenum != 1)
  {
    fprintf(stderr, "error: expected 1 input file, got %d\n", filenum);
    print_usage(stderr);
    return 1;
  }
  const char *infile = argv[optind];

  // Add all taxa to a hashtable
  FILE *taxnodes = myfopen(nodesfile, "r");
  Hashmap *nodes = hashmap_new(500000);
  char buffer[1024];
  while( (fgets(buffer, 1024, taxnodes)) != NULL )
  {
    char *nid = strtok(buffer, "|");
    char *pid = strtok(NULL, "|");
    char *rank = strtok(NULL, "|");

    nid = trim_str(nid);
    pid = trim_str(pid);
    rank = trim_str(rank);

    Taxon *t = taxon_new(nid, pid, rank);
    int store_success = hashmap_add(nodes, nid, t);
    if(store_success == 0)
    {
      fprintf(stderr, "error storing %p (id=%s) to hashmap\n", t, t->id);
      return 1;
    }
  }
  fclose(taxnodes);

  // Resolve parent/child relationships
  char **tids = hashmap_keys(nodes);
  int i;
  for(i = 0; i < hashmap_size(nodes); i++)
  {
    Taxon *t = hashmap_get(nodes, tids[i]);
    if(t == NULL)
    {
      fprintf(stderr, "error fetching taxon with tid '%s'\n", tids[i]);
      return 1;
    }

    Taxon *parent = hashmap_get(nodes, t->parentid);
    if(parent == NULL)
    {
      fprintf(stderr, "error fetching parent with id '%s'\n", t->parentid);
      return 1;
    }
    t->parent = parent;
  }
  free(tids);

  // Add names to taxa
  FILE *taxnames = myfopen(namesfile, "r");
  while( (fgets(buffer, 1024, taxnames)) != NULL )
  {
    char *tid = strtok(buffer, "|");
    char *name = strtok(NULL, "|");
    char *uniqname = strtok(NULL, "|");
    char *class = strtok(NULL, "|");

    tid = trim_str(tid);
    name = trim_str(name);
    uniqname = trim_str(uniqname);
    class = trim_str(class);

    Taxon *t = hashmap_get(nodes, tid);
    if(t == NULL)
    {
      fprintf(stderr, "error fetching node '%s|%s' for naming\n", tid, name);
      return 1;
    }
    if(strcmp(class, "scientific name") == 0)
      t->name = strdup(name);
  }
  fclose(taxnames);

  // Determine classification of each taxon at the specified rank
  FILE *instream = myfopen(infile, "r");
  while( (fgets(buffer, 1024, instream)) != NULL )
  {
    char *taxonid = strtok(buffer, delims);
    Taxon *taxon = hashmap_get(nodes, taxonid);
    if(taxon == NULL)
    {
      fprintf(stderr, "error fetching taxon by id '%s'\n", taxonid);
      return 1;
    }

    Taxon *target = get_parent_by_rank(taxon, targetrank);
    if(target == NULL)
    {
      fprintf(stderr, "error resolving %s for species %s\n", targetrank, taxonid);
      return 1;
    }
    else
    {
      if(filterrank != NULL && filtertaxon != NULL)
      {
        Taxon *filtertarget = get_parent_by_rank(taxon, filterrank);
        int matchfilter = strcmp(filtertarget->name, filtertaxon) == 0;
        if(!(matchfilter ^ complement))
          continue;
      }

      printf("%s,%s,%s\n", taxonid, target->id, target->name);
    }
  }
  fclose(instream);

  hashmap_delete(nodes, (void (*)(void *))taxon_delete);
  return 0;
}
