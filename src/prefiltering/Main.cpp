#include <iostream>
#include <unistd.h>
#include <string>
#include <signal.h>
#include <execinfo.h>
#include "Prefiltering.h"
#include "CommandDeclarations.h"


#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif


void mmseqs_debug_catch_signal(int sig_num)
{
  if(sig_num == SIGILL)
  {
    fprintf(stderr, "Your CPU does not support all the latest features that this version of mmseqs makes use of!\n"
                    "Please run on a newer machine.");
    EXIT(sig_num);
  }
  else
  {
    fprintf (stderr, "\n\n-------------------------8<-----------------------\nExiting on error!\n");
    fprintf (stderr, "Signal %d received\n",sig_num);
    perror("ERROR (can be bogus)");

    fprintf(stderr, "Backtrace:");
    void *buffer[30];
    int nptrs = backtrace(buffer, 30);
    backtrace_symbols_fd(buffer, nptrs, 2);
    fprintf (stderr, "------------------------->8-----------------------\n\n"
        "Send the binary program that caused this error and the coredump (ls core.*).\n"
        "Or send the backtrace:"
        "\n$ gdb -ex=bt --batch PROGRAMM_NAME CORE_FILE\n"
        "If there is no core file, enable coredumps in your shell and run again:\n"
        "$ ulimit -c unlimited\n\n");
  }

  EXIT(1);
}

void mmseqs_cuticle_init()
{
  struct sigaction handler;
  handler.sa_handler = mmseqs_debug_catch_signal;
  sigemptyset(&handler.sa_mask);
  handler.sa_flags = 0;

  sigaction(SIGFPE, &handler, NULL);
  sigaction(SIGSEGV, &handler, NULL);
  sigaction(SIGBUS, &handler, NULL);
  sigaction(SIGABRT, &handler, NULL);
}

void printUsagePrefiltering(){

    std::string usage("\nCalculates k-mer similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de) & Martin Steinegger (Martin.Steinegger@campus.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_pref <queryDB> <targetDB> <outDB> [opts]\n"
            "-s              \t[float]\tSensitivity in the range [1:9] (default=4).\n"
            "-k              \t[int]\tk-mer size in the range [4:7] (default=6).\n"
            "-cpu              \t[int]\tNumber of cores used for the computation (default=all cores).\n"
            "--alph-size     \t[int]\tAmino acid alphabet size (default=21).\n"
            "--z-score       \t[float]\tZ-score threshold (default: 50.0).\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n"
            "--profile       \t\tHMM Profile input.\n"
            "--nucl          \t\tNucleotide sequences input.\n"
            "--search-mode   \t[int]\tSearch mode loc: 1 glob: 2 (default=1).\n"
            "--no-comp-bias-corr \t\tSwitch off local amino acid composition bias correction.\n"
            "--no-spaced-kmer \t\tSwitch off spaced kmers (consecutive pattern).\n"
            "--split         \t[int]\tSplits target databases in n equal distrbuted junks (default=1)\n"
            "--threads       \t[int]\tNumber of threads used to compute. (Default=all cpus)\n"
            "--max-seqs      \t[int]\tMaximum result sequences per query (default=300).\n"
            "--skip          \t[int]\tNumber of skipped k-mers during the index table generation.\n"
            "--sub-mat       \t[file]\tAmino acid substitution matrix file.\n"
            "-v              \t[int]\tVerbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO (default=3).\n");
    Debug(Debug::INFO) << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase,
               std::string* ffindexTargetDBBase, std::string* ffindexOutDBBase,
               std::string* scoringMatrixFile, float* sens, int* kmerSize,
               int* alphabetSize, float* zscoreThr, size_t* maxSeqLen,
               int* querySeqType, int * targetSeqType, bool* localSearch,
               size_t* maxResListLen, bool* compBiasCorrection,
               bool* spacedKmer, int* split, int* threads, int* skip, int* verbosity){

    if (argc < 4){
        printUsagePrefiltering();
        EXIT(EXIT_FAILURE);
    }

    ffindexQueryDBBase->assign(argv[1]);
    ffindexTargetDBBase->assign(argv[2]);
    ffindexOutDBBase->assign(argv[3]);

    int zscoreSet = 0;

    int i = 4;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
            if (*querySeqType == Sequence::NUCLEOTIDES){
                Debug(Debug::ERROR) << "No scoring matrix is allowed for nucleotide sequences.\n";
                EXIT(EXIT_FAILURE);
            }
            if (++i < argc){
                scoringMatrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        } 
        else if (strcmp(argv[i], "--sub-mat") == 0){
            if (*querySeqType == Sequence::NUCLEOTIDES){
                Debug(Debug::ERROR) << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            if (++i < argc){
                scoringMatrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *sens = atof(argv[i]);
                if (*sens < 1.0 || *sens > 9.0){
                    Debug(Debug::ERROR) << "Please choose sensitivity in the range [1:9].\n";
                    exit(EXIT_FAILURE);
                }
                // adapt z-score threshold to the sensitivity setting
                // user defined threshold overwrites the automatic setting
                if (zscoreSet == 0){
                    if (1.0 <= *sens && *sens < 2.0)
                        *zscoreThr = 500.0;
                    else if (2.0 <= *sens && *sens < 3.0)
                        *zscoreThr = 300.0;
                    else if (3.0 <= *sens && *sens < 4.0)
                        *zscoreThr = 100.0;
                    else if (4.0 <= *sens && *sens < 5.0)
                        *zscoreThr = 50.0;
                    else if (5.0 <= *sens && *sens < 6.0)
                        *zscoreThr = 40.0;
                    else if (6.0 <= *sens && *sens < 7.0)
                        *zscoreThr = 30.0;
                    else if (7.0 <= *sens && *sens < 8.0)
                        *zscoreThr = 20.0;
                    else if (8.0 <= *sens && *sens <= 9.0)
                        *zscoreThr = 10.0;
                }
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-k") == 0){
            if (++i < argc){
                *kmerSize = atoi(argv[i]);
                if (*kmerSize < 4 || *kmerSize > 7){
                    Debug(Debug::ERROR) << "Please choose k in the range [4:7].\n";
                    EXIT(EXIT_FAILURE);
                }
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-a") == 0){
            if (++i < argc){
                *alphabetSize = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--alph-size") == 0){
            if (++i < argc){
                *alphabetSize = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--z-score-thr") == 0){
            if (++i < argc){
                *zscoreThr = atof(argv[i]);
                zscoreSet = 1;
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided" << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--z-score") == 0){
            if (++i < argc){
                *zscoreThr = atof(argv[i]);
                zscoreSet = 1;
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided" << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--max-seq-len") == 0){
            if (++i < argc){
                *maxSeqLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--nucl") == 0){
            if (strcmp(scoringMatrixFile->c_str(), "") != 0){
                Debug(Debug::ERROR) << "No scoring matrix is allowed for nucleotide sequences.\n";
                EXIT(EXIT_FAILURE);
            }                                                           
            *querySeqType = Sequence::NUCLEOTIDES;
            *targetSeqType = Sequence::NUCLEOTIDES;
            i++;
        }
        else if (strcmp(argv[i], "--profile") == 0){
            *querySeqType = Sequence::HMM_PROFILE;
            *targetSeqType = Sequence::AMINO_ACIDS;

            i++;
        }
        else if (strcmp(argv[i], "--max-seqs") == 0){
            if (++i < argc){
                *maxResListLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--no-comp-bias-corr") == 0){
            *compBiasCorrection = false;
            i++;
        }
        else if (strcmp(argv[i], "--no-spaced-kmer") == 0){
            *spacedKmer = false;
            i++;
        }
        else if (strcmp(argv[i], "--search-mode") == 0){
            if (++i < argc){
                switch (atoi(argv[i])) {
                    default:
                    case 1:
                        *localSearch = true;
                        break;
                    case 2:
                        *localSearch = false;
                        break;
                }
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--skip") == 0){
            if (++i < argc){
                *skip = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-cpu") == 0){
            if (++i < argc){
                *threads = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-v") == 0){
            if (++i < argc){
                *verbosity = atoi(argv[i]);
                if (*verbosity < 0 || *verbosity > 3){
                    Debug(Debug::ERROR) << "Wrong value for verbosity, please choose one in the range [0:3].\n";
                    EXIT(1);
                }
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--split") == 0){
            if (++i < argc){
                *split = atoi(argv[i]);
                i++;
            }
            else {
                printUsagePrefiltering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        else {
            printUsagePrefiltering();
            Debug(Debug::ERROR) << "Wrong argument: " << argv[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}



int prefilter(int argc, const char **argv)
{
#ifdef HAVE_MPI 
    int mpi_error,mpi_rank,mpi_num_procs;
    mpi_error = MPI_Init(&argc, &argv);
    mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    mpi_error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
    Debug(Debug::WARNING) << "MPI Init...\n";
    Debug(Debug::WARNING) << "Rank: " << mpi_rank << " Size: " << mpi_num_procs << "\n";


#endif
    mmseqs_cuticle_init();

    int verbosity = Debug::INFO;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    int kmerSize =  6;
    int alphabetSize = 21;
    size_t maxSeqLen = 50000;
    size_t maxResListLen = 300;
    float sensitivity = 4.0f;
    int split = 1;
    int skip = 0;
    int querySeqType  = Sequence::AMINO_ACIDS;
    int targetSeqType = Sequence::AMINO_ACIDS;

    int threads = 1;
#ifdef OPENMP
    threads = Util::omp_thread_count();
#endif
    bool compBiasCorrection = true;
    bool spacedKmer = true;
    bool localSearch = true;


    float zscoreThr = 50.0f;

    std::string queryDB = "";
    std::string targetDB = "";
    std::string outDB = "";
    // get the path of the scoring matrix
    char* mmdir = getenv ("MMDIR");
    if (mmdir == 0){
        std::cerr << "Please set the environment variable $MMDIR to your MMSEQS installation directory.\n";
        exit(1);
    }
    std::string scoringMatrixFile(mmdir);
    scoringMatrixFile = scoringMatrixFile + "/data/blosum62.out";

    // print command line
    Debug(Debug::WARNING) << "Program call:\n";
    for (int i = 0; i < argc; i++)
        Debug(Debug::WARNING) << argv[i] << " ";
    Debug(Debug::WARNING) << "\n\n";

    parseArgs(argc, argv, &queryDB, &targetDB, &outDB, &scoringMatrixFile,
                          &sensitivity, &kmerSize, &alphabetSize, &zscoreThr,
                          &maxSeqLen, &querySeqType, &targetSeqType, &localSearch,
                          &maxResListLen, &compBiasCorrection,
                          &spacedKmer, &split, &threads, &skip, &verbosity);

#ifdef OPENMP
    omp_set_num_threads(threads);
#endif

    Debug::setDebugLevel(verbosity);

    if (querySeqType == Sequence::NUCLEOTIDES)
        alphabetSize = 5;

    Debug(Debug::WARNING) << "k-mer size: " << kmerSize << "\n";
    Debug(Debug::WARNING) << "Alphabet size: " << alphabetSize << "\n";
    Debug(Debug::WARNING) << "Sensitivity: " << sensitivity << "\n";
    Debug(Debug::WARNING) << "Z-score threshold: " << zscoreThr << "\n";
    if (compBiasCorrection)
        Debug(Debug::WARNING) << "Compositional bias correction is switched on.\n";
    else
        Debug(Debug::WARNING) << "Compositional bias correction is switched off.\n";
    if (spacedKmer)
        Debug(Debug::WARNING) << "Spaced kmers is switched on.\n";
    else
        Debug(Debug::WARNING) << "Spaced kmers is switched off.\n";
    Debug(Debug::WARNING) << "\n";

    std::string queryDBIndex = queryDB + ".index";
    std::string targetDBIndex = targetDB + ".index";
    std::string outDBIndex = outDB + ".index";

    Debug(Debug::WARNING) << "Initialising data structures...\n";
    Prefiltering* pref = new Prefiltering(queryDB, queryDBIndex, targetDB, targetDBIndex,
                                          outDB, outDBIndex, scoringMatrixFile, sensitivity,
                                          kmerSize, spacedKmer, maxResListLen, alphabetSize, zscoreThr,
                                          maxSeqLen, querySeqType, targetSeqType, compBiasCorrection,
                                          split, skip);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    gettimeofday(&start, NULL);
#ifdef HAVE_MPI
    pref->run(mpi_rank, mpi_num_procs);
#else
    pref->run();
#endif
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "\nTime for prefiltering run: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    delete pref;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
