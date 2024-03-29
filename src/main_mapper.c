/*
    intmap -p file
    

*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
//#include <time.h>

#include "version.h"
#include "main_mapper.h"
#include "genfa.h"
#include "genfq.h"
#include "runall.h"
#include "tech_map.h"
#include "extract.h"

#define TRUE  1
#define FALSE 0


/*
 * 
 * */

#define USAGE_ITR 0
#define USAGE_CMD 1
#define USAGE_OPT 2


extern void usage_runall( int opt )
{
    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "runall : Proceed genfqse/genfqpe/analyze\n" ); 
    }

    if( opt == USAGE_OPT )
    {
        fprintf( stderr, "    -l   Length [1-200]\n" ); 
        fprintf( stderr, "    -q   Quality  value. Refer http://en.wikipedia.org/wiki/FASTQ_format\n" ); 
        fprintf( stderr, "    -s   Substitution rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -i   Insertion rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -d   Deletion rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -o   library distance bp [100-250] between end of start and start of end\n" ); 
        fprintf( stderr, "         If -o option is specified, pair-end aligment will be forced unless it is 0(zero)\n" ); 
        fprintf( stderr, "         If -o option is zero, single-end reads will be generated\n" ); 
        fprintf( stderr, "    -t   Threshold for map quality to analyze sam file and generatemap result.\n" ); 
        fprintf( stderr, "    -f   FASTA file name to be generated, e.g. ref.fa\n" ); 
        fprintf( stderr, "    -b   The number of bases per line in FASTA file before it is preprocessed\n" ); 
        fprintf( stderr, "    -x   index file path\n" ); 
        fprintf( stderr, "    -p   use if you use Hadoop and should specify the directory that includes bwa and samtools,\n" ); 
        fprintf( stderr, "         and can be accessed by hadoop.\n" ); 
        fprintf( stderr, "\n" );
    }
}

extern void usage_genfa( int opt )
{
    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "genfa : Generate FastA file from preprocessed file\n" ); 
    }
    if( opt == USAGE_OPT )
    {
        fprintf( stderr, "    -f   FASTA file name to be generated, e.g. ref.fa\n" ); 
        fprintf( stderr, "    -b   The number of bases in a line of FASTA file\n" ); 
        fprintf( stderr, "\n" );
    }
}

extern void usage_genfqse( int opt )
{
    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "genfqse : Generate FastQ file for single ends\n" ); 
    }

    if( opt == USAGE_OPT )
    {
        fprintf( stderr, "    -l   Length [1-200]\n" ); 
        fprintf( stderr, "    -q   Quality  value. Refer http://en.wikipedia.org/wiki/FASTQ_format\n" ); 
        fprintf( stderr, "    -s   Substitution rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -i   Insertion rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -d   Deletion rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "\n" );
    }
}


extern void usage_genfqpe( int opt )
{
    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "genfqpe : Generate FastQ file for mate pair\n" ); 
    }
    if( opt == USAGE_OPT )
    {
        fprintf( stderr, "    -l   Length [1-200]\n" ); 
        fprintf( stderr, "    -q   Quality  value. Refer http://en.wikipedia.org/wiki/FASTQ_format\n" ); 
        fprintf( stderr, "    -s   Substitution rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -i   Insertion rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -d   Deletion rate. [0.0000 - 1.0000] 1 means 100%%, 0.25 means 25%%\n" ); 
        fprintf( stderr, "    -o   library distance bp [100-250] between end of start and start of end\n" ); 
        fprintf( stderr, "         If -o option is specified, pair-end aligment will be forced unless it is 0(zero)\n" ); 
        fprintf( stderr, "         If -o option is zero, single-end reads will be generated\n" ); 
        fprintf( stderr, "    -b   The number of bases per line in FASTA file before it is preprocessed\n" ); 
        fprintf( stderr, "\n" );
    }
}

extern void usage_extract( int opt )
{
    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "extract : extract information from bam file then make an intermediate file\n" ); 
    }
    if( opt == USAGE_OPT )
    {
        fprintf( stderr, "    -l   Length [1-200]\n" ); 
        fprintf( stderr, "    -f   FASTA file name to be generated, e.g. ref.fa\n" ); 
        fprintf( stderr, "    -m   BAM file name\n" ); 
        fprintf( stderr, "\n" );
    }
}


extern void usage_tech( int opt )
{
    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "tech : automatically set variables by sequencing technology\n" ); 
    }
    if( opt == USAGE_OPT )
    {
        fprintf( stderr, "    --illumina | --solid | --roche | pacbio | pacbio_ec\n" ); 
        fprintf( stderr, "    -b   The number of bases per line in FASTA file before it is preprocessed\n" ); 
        fprintf( stderr, "    -x   index file path\n" ); 
        fprintf( stderr, "\n" );
    }
}


extern void usage( int opt )
{
    if( opt == USAGE_ITR || opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "\nProgram : GMA (Genome Mappability Analyzer) MAPPER\n" );
        fprintf( stderr, "Version : %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD );
        fprintf( stderr, "Usage : gma [COMMAND] [OPTIONS] [FASTA]\n\n");
    }

    if( opt == USAGE_CMD || opt == USAGE_OPT )
    {
        fprintf( stderr, "[COMMAND] runall | genfa | genfqse | genfqpe | extract | tech \n" );
        fprintf( stderr, "\n" );
    }
 
    if( opt == USAGE_OPT )
    {
        usage_runall( opt );
        usage_genfa( opt );
        usage_genfqse( opt );
        usage_genfqpe( opt );
        usage_extract( opt );
        usage_tech( opt );
    }
}


/*
 *
 */

#define ARGC_MAX 50

int main( int argc, char *argv[] )
{
    
    char *ARGV[ARGC_MAX];
    int i = 0;
    int ARGC = 0;
    int ret = EXIT_SUCCESS;

    if ( argc == 1 )    // HADOOP
    {
        ARGV[0] = getenv("mapred_job_arg0"); // # of args, i.e., argc
        ARGV[1] = getenv("mapred_job_arg1"); // e.g. -l
        ARGV[2] = getenv("mapred_job_arg2"); // e.g. -l
        ARGV[3] = getenv("mapred_job_arg3"); // e.g. 50
        ARGV[4] = getenv("mapred_job_arg4"); // e.g. -e
        ARGV[5] = getenv("mapred_job_arg5"); // e.g. 0.02
        ARGV[6] = getenv("mapred_job_arg6"); // e.g. -q
        ARGV[7] = getenv("mapred_job_arg7"); // e.g. A
        ARGV[8] = getenv("mapred_job_arg8"); // e.g. -r
        ARGV[9] = getenv("mapred_job_arg9"); // e.g. 200
        ARGV[10] = getenv("mapred_job_arg10"); // e.g. -t
        ARGV[11] = getenv("mapred_job_arg11"); // e.g. 20
        ARGV[12] = getenv("mapred_job_arg12"); // e.g. -f
        ARGV[13] = getenv("mapred_job_arg13"); // e.g. ref.fa
        ARGV[14] = getenv("mapred_job_arg14"); // e.g. ref.fa
        ARGV[15] = getenv("mapred_job_arg15"); // e.g. ref.fa
        ARGV[16] = getenv("mapred_job_arg16"); // e.g. ref.fa
        ARGV[17] = getenv("mapred_job_arg17"); // e.g. ref.fa
        ARGV[18] = getenv("mapred_job_arg18"); // e.g. ref.fa
        ARGV[19] = getenv("mapred_job_arg19"); // e.g. ref.fa
        ARGV[20] = getenv("mapred_job_arg20"); // e.g. ref.fa
        ARGV[21] = getenv("mapred_job_arg21"); // e.g. ref.fa
        ARGV[22] = getenv("mapred_job_arg22"); // e.g. ref.fa
        ARGV[23] = getenv("mapred_job_arg23"); // e.g. ref.fa
       
       
        if( ARGV[0] == NULL )
        {
            usage( USAGE_CMD );
        }
        else
        { 
            fprintf( stderr, "=============================================\n");
            fprintf( stderr, "[HADOOP] # of args : %s\n", ARGV[0] );
            fprintf( stderr, "[HADOOP] command : %s\n", ARGV[1] );
    
            ARGC = atoi( ARGV[0] );
    
            for ( i = 0; i < ARGC; i++ )
            {
                fprintf( stderr, "%d. %s\n", i , ARGV[i] );
            }
            fprintf( stderr, "=============================================\n");
    
            // For hadoop, only "runall" is allowed
    
            if( strcmp( ARGV[1], "runall") == 0 ) ret = runall( ARGC-1, ARGV+1 );
            else if( strcmp( ARGV[1], "genfa") == 0 ) genfa( ARGC-1, ARGV+1 );
            else if( strcmp( ARGV[1], "genfqse") == 0 ) genfq( ARGC-1, ARGV+1 );
            else if( strcmp( ARGV[1], "genfqpe") == 0 ) genfq( ARGC-1, ARGV+1 );
            else if( strcmp( ARGV[1], "extract") == 0 ) extract( ARGC-1, ARGV+1 );
            else if( strcmp( ARGV[1], "tech") == 0 ) runall_by_tech( ARGC-1, ARGV+1 );
            else usage( USAGE_OPT );
        }
    }
    else if( argc == 2 ) 
    {
        if( strcmp( "-h", argv[1] ) == 0)
        {
             usage( USAGE_OPT );
        }
        else if ( strcmp( "-v", argv[1] ) == 0)  
        {
            fprintf( stderr, "%d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD );
        }
        else
        {
            usage( USAGE_CMD );

            if( strcmp( argv[1], "runall") == 0 ) usage_runall( USAGE_OPT );
            else if( strcmp( argv[1], "genfa") == 0 ) usage_genfa( USAGE_OPT );
            else if( strcmp( argv[1], "genfqse") == 0 ) usage_genfqse( USAGE_OPT );
            else if( strcmp( argv[1], "genfqpe") == 0 ) usage_genfqpe( USAGE_OPT );
            else if( strcmp( argv[1], "extract") == 0 ) usage_extract( USAGE_OPT );
            else if( strcmp( argv[1], "tech") == 0 ) usage_tech( USAGE_OPT );
            else usage( USAGE_OPT );
        }
    }
    else
    { 
        fprintf( stderr, "=============================================\n");
        fprintf( stderr, "[LOCAL] # of args : %d\n", argc );
        fprintf( stderr, "[LOCAL] command : %s\n", argv[0] );
        
        for ( i = 1; i < argc; i++ )
        {
            fprintf( stderr, "%d. %s\n", i, argv[i] );
        }

        fprintf( stderr, "=============================================\n");

        if( strcmp( argv[1], "runall") == 0 ) runall( argc-1, argv+1 );
        else if( strcmp( argv[1], "genfa") == 0 ) genfa( argc-1, argv+1 );
        else if( strcmp( argv[1], "genfqse") == 0 ) genfq( argc-1, argv+1 );
        else if( strcmp( argv[1], "genfqpe") == 0 ) genfq( argc-1, argv+1 );
        else if( strcmp( argv[1], "extract") == 0 ) extract( argc-1, argv+1 );
        else if( strcmp( argv[1], "tech") == 0 ) runall_by_tech( argc-1, argv+1 );
        else usage( USAGE_OPT );
        
    }

    return ret;
                          //EXIT_SUCCESS;  // must return 0(=EXIT_SUCCESS)
                          // non-zero means errors, EXIT_FAILURE
}


