#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

/*
 * The majority of the sequential code in this project is replicated across..-
 * -.. all versions of this algorithm - MPI, Pthreads and OpenMP.
 */


struct SharedData {
    char*** patterns;
    char** sequences;
    int** foundMatches;
    size_t* patternLengths;
    int* patternCount;
    int* matchCounter;
};

struct ThreadData {
    int start;
    int end;
    int extraThreads;
    int sequenceCount;
    pthread_barrier_t syncBarrier;
    struct SharedData sharedData;
};

struct InnerThreadData {
    char** currentPatterns;
    char* commonStart;
    size_t* currentSize;
    int start;
    int end;
    int portion;
    int currIndex;
    struct SharedData sharedData;
    pthread_barrier_t syncBarrier;
    pthread_mutex_t resizeMutex;
    pthread_mutex_t addMatchMutex;
};


int getSequenceData(FILE* file, int* sequenceCount, int* maxSequenceLength) {
    fscanf(file, "%d", sequenceCount);
    fscanf(file, "%d", maxSequenceLength);

    if (!sequenceCount || sequenceCount < 0 ||
        !maxSequenceLength || maxSequenceLength < 0) {
        printf("Invalid input.\n");
        return 0;
    }
    return 1;
}

int malloc_pattern_arrays(char**** patterns, size_t** patternLengths,
                          int elementCount) {
    *patterns = (char***) malloc(elementCount * sizeof(char**));
    if (*patterns) {
        *patternLengths = (size_t *) malloc(elementCount * sizeof(size_t));
        if (!*patternLengths) {
            free(*patterns);
            return 0;
        }
        return 1;
    }
    return 0;
}

void free_inner_2DArray(void** array, size_t length) {
    for (int i = 0; i < length; i++) {
        free(array[i]);
    }
}

void free_full_2DArray(void** array, int length) {
    free_inner_2DArray(array, length);
    free(array);
}


int malloc_pattern_array(char*** array1, size_t elementCount, int repeat) {
    // The function parameters are specific for less complexity as this is the only purpose they are used for.
    for (int i = 0; i < repeat; i++) {
        array1[i] = (char**) malloc(elementCount * sizeof(char *));
        if (!array1[i]) {
            free_full_2DArray((void **) array1, i);
            return 0;
        }
    }
    return 1;
}

void free_pattern_arrays(char*** patterns, size_t* patternLengths,
                         int repeats, int size, int* patternCount) {
    for (int i = 0; i < repeats; i++) {
        free_inner_2DArray((void**) patterns[i], patternCount[i]);
    }
    free_full_2DArray((void**) patterns, size);
    free(patternLengths);
}

FILE* getFile(char* fileName) {
    FILE* file = fopen(fileName, "r");
    return file;
}

int distinguishPatterns(FILE* file, char*** patterns, size_t* patternLengths,
                        int sequenceCount, size_t estimatedPatternCount,
                        int* patternCount) {
    char patternBuffer[255]; // Max length of a pattern
    int activeBracket = 0;
    int k;
    int currentVariations;
    size_t currentLength;
    char currentLetter;
    for (int i = 0; i < sequenceCount; i++){
        fscanf(file, "%s", patternBuffer);
        currentLength = strlen(patternBuffer);
        char withoutVariation[currentLength]; // Pattern before variable input
        currentVariations = 0;
        k = 0;
        for (int j = 0; j < currentLength; j++) {
            currentLetter = patternBuffer[j];
            if (currentLetter == '[') {
                // Embedded brackets or more than one variation in the same pattern
                if (activeBracket || currentVariations) {
                    printf("Pattern %d is invalid.\nEmbedded brackets are not"
                           " allowed and only one variation is permitted"
                           " within each pattern.\n", i);
                    free_pattern_arrays(patterns, patternLengths, i,
                                        sequenceCount, patternCount);
                    return 0;
                }
                withoutVariation[k] = '\0';
                activeBracket++;
                continue;
            }
            if (currentLetter == ']') {
                // No corresponding open bracket, or empty bracket
                if (!activeBracket || !currentVariations) {
                    printf("Pattern %d is invalid.\nClosing brackets "
                           "must have corresponding opening brackets and "
                           "brackets may not be empty.\n", i);
                    free_pattern_arrays(patterns, patternLengths, i,
                                        sequenceCount, patternCount);
                    return 0;
                }
                // Variation within current pattern already exists - limited to 1.
                activeBracket--;
                k = 0;
                continue;
            }
            if (activeBracket) {
                if (currentVariations == estimatedPatternCount) {
                    patterns[i] = (char**) realloc(
                            patterns[i],currentVariations*2 * sizeof(char*));

                    if (!patterns[i] || patterns[i][currentVariations]) {
                        free_pattern_arrays(patterns, patternLengths, i,
                                            sequenceCount, patternCount);
                        return 0;
                    }
                }
                // currentLength is either the whole length or includes 2 brackets.
                patterns[i][currentVariations] = (char*)
                        malloc((currentLength-2+1) * sizeof(char));
                if (!patterns[i][currentVariations]) {
                    free_pattern_arrays(patterns, patternLengths, i,
                                        sequenceCount, patternCount);
                    return 0;
                }
                // Store one variation of the pattern
                strncpy(patterns[i][currentVariations], withoutVariation, k);
                patterns[i][currentVariations][k] = '\0';
                strncat(patterns[i][currentVariations], &currentLetter, 1);
                patternCount[i]++;
                currentVariations++;
            }
            else {
                withoutVariation[k] = currentLetter;
                k++;
            }
        }
        withoutVariation[k] = '\0';
        if (!currentVariations) {
            patterns[i] = (char**) realloc(patterns[i], 1 * sizeof(char*));
            if (patterns[i]) {
                patterns[i][0] = (char*) malloc((k+1) * sizeof(char));
            }
            if (!patterns[i] || !patterns[i][0]) {
                free_pattern_arrays(patterns, patternLengths, i,
                                    sequenceCount, patternCount);
                return 0;
            }
            strncpy(patterns[i][0], withoutVariation, k+1);
            patternLengths[i] = k;
            patternCount[i] = 1;
            continue;
        }
        currentLength = 0;
        for (int j = 0; j < currentVariations; j++) {
            strncat(patterns[i][j], withoutVariation, k+1);
            if (!currentLength) { // All patterns variations are of the same length
                currentLength = strlen(patterns[i][j]);
            }
            patternLengths[i] = currentLength;
            patterns[i][j] = (char*) realloc(patterns[i][j],
                                             (currentLength + 1) * sizeof(char));
            if (!patterns[i][j]) {
                free_pattern_arrays(patterns, patternLengths, i,
                                    sequenceCount, patternCount);
                return 0;
            }
        }
        patternCount[i] = currentVariations;
    }
    return 1;
}

int parseSequences(FILE* file, char** sequences, int sequenceCount, int maxSequenceLength) {
    char seqBuffer[maxSequenceLength+1];
    size_t currentLength;
    for (int i = 0; i < sequenceCount; i++) {
        fscanf(file, "%s", seqBuffer);
        currentLength = strlen(seqBuffer);
        sequences[i] = (char*) malloc(currentLength * sizeof(char));
        if (!sequences[i]) {
            free_full_2DArray((void**) sequences, i); // Free all up to i
            return 0;
        }
        seqBuffer[currentLength] = '\0';
        strncpy(sequences[i], seqBuffer, currentLength+1);
    }
    return 1;
}

void findCommonStart(char* commonStart, char** patterns, size_t patternLength,
                     int index, const int* patternCount) {
    for (int j = 0; j < patternLength; j++) {
        for (int k = 1; k < patternCount[index]; k++) {
            if (patterns[k][j] != patterns[k - 1][j]) {
                commonStart[j] = '\0';
                j = patternLength; // End outer loop - factor whole loop into function later and return instead
                break;
            }
            commonStart[j] = patterns[0][j]; // index 0 is allowed because all instances of the first letter are the same per the above loop
        }
    }
}

int patternMatchCommonStart(int patternCount, size_t patternLength,
                            const char* matchedString, char** patterns) {
    for (int j = 0; j < patternCount; j++) {
        for (int k = 0; k < patternLength; k++) {
            if (matchedString[k] != patterns[j][k]) {
                break;
            }
            if (k == patternLength-1) return 1; // Found
        }
    }
    return 0;
}

int findMatches(char* sequence, char** currentPatterns, char* commonStart,
                int** foundMatches, int* currMatchCounter, size_t* currentSize,
                int patternCount, int currIndex, size_t patternLength, int offset,
                pthread_mutex_t* resizeMutex, pthread_mutex_t* addMatchMutex,
                int parallel) {
    char* temp = sequence;
    int found, index;
    while (temp[0] != '\0') { // End of sequence
        char* commonStartMatch = strstr(temp, commonStart); // Shortcut - initially only look for common start
        if (!commonStartMatch) break; // Not found
        found = patternMatchCommonStart(patternCount, patternLength, commonStartMatch,
                                        currentPatterns);
        if (found) {
            if (parallel) {
                pthread_mutex_lock(addMatchMutex);
                index = (int) (commonStartMatch - sequence) + offset;
                foundMatches[currIndex][*currMatchCounter] = index;
                (*currMatchCounter)++;
                pthread_mutex_unlock(addMatchMutex);
            }
            else {
                foundMatches[currIndex][*currMatchCounter] = (int) (commonStartMatch - sequence);
                (*currMatchCounter)++;
            }
            temp = &(temp[commonStartMatch - temp + 1]);
        }
        else temp = &temp[1]; // Skip one character to remove the already-found common start
        if (*currMatchCounter >= *currentSize) { // Enlarge array
            if (parallel) {
                pthread_mutex_lock(resizeMutex);
                if (*currMatchCounter < *currentSize) {
                    pthread_mutex_unlock(resizeMutex);
                    continue;
                }
                (*currentSize) *= 2;
                foundMatches[currIndex] = (int *)
                        realloc(foundMatches[currIndex], (*currentSize) * sizeof(int));
                pthread_mutex_unlock(resizeMutex);
            }
            else {
                (*currentSize) *= 2;
                foundMatches[currIndex] = (int*)
                        realloc(foundMatches[currIndex], (*currentSize) * sizeof(int));
            }
            if (!foundMatches[currIndex]) {
                return 0;
            }
        }
    }
    return 1;
}

int parallel_manageMatches(struct InnerThreadData data) {
    int start = data.start;
    int end = data.end;
    int portion = data.portion;
    pthread_barrier_wait(&data.syncBarrier);

    int currIndex = data.currIndex;
    size_t patternLength = data.sharedData.patternLengths[currIndex];
    char *sequence = data.sharedData.sequences[currIndex];
    int **foundMatches = data.sharedData.foundMatches;
    int *matchCounter = data.sharedData.matchCounter;

    // Start earlier in case split is at a matching case
    if (start) { // Not the first chunk
        start -= (int) patternLength - 1;
        portion += (int) patternLength - 1;
    }
    char *localSequence = (char *) malloc(portion + 1 * sizeof(char));
    if (!localSequence) return 0;

    memcpy(localSequence, &(sequence[start]), portion * sizeof(char));
    localSequence[portion] = '\0';

    if (!findMatches(localSequence, data.currentPatterns, data.commonStart,
                     foundMatches, &matchCounter[currIndex], data.currentSize,
                     data.sharedData.patternCount[currIndex], currIndex,
                     patternLength, start, &data.resizeMutex,
                     &data.addMatchMutex, 1)) {
        free(localSequence);
        return 0;
    }
    free(localSequence);
    return 1;
}

int manageMatches(struct ThreadData data) {
    int start = data.start;
    int end = data.end;
    pthread_barrier_wait(&data.syncBarrier);
    char*** patterns = data.sharedData.patterns;
    char** sequences = data.sharedData.sequences;
    int** foundMatches = data.sharedData.foundMatches;
    size_t* patternLengths = data.sharedData.patternLengths;
    int* patternCount = data.sharedData.patternCount;
    int* matchCounter = data.sharedData.matchCounter;
    int threadsPerThread = data.extraThreads;

    size_t currentSizes[data.sequenceCount];

    pthread_t* threads;
    struct InnerThreadData threadData;
    if (threadsPerThread > 1) {
        threadData.sharedData = data.sharedData;

        threads = (pthread_t*) malloc(threadsPerThread * sizeof(pthread_t));
        if (!threads) return 0;

        if (pthread_barrier_init(&threadData.syncBarrier, NULL, 2)) {
            free(threads);
            return 0;
        }

        if (pthread_mutex_init(&threadData.resizeMutex, NULL)) {
            free(threads);
            pthread_barrier_destroy(&threadData.syncBarrier);
            return 0;
        }

        if (pthread_mutex_init(&threadData.addMatchMutex, NULL)) {
            free(threads);
            pthread_barrier_destroy(&threadData.syncBarrier);
            pthread_mutex_destroy(&threadData.resizeMutex);
            return 0;
        }
    }

    char** currentPatterns;
    size_t patternLength, currentSize;
    for (int i = start; i < end; i++) {
        currentSizes[i] = 3; // Estimate 3 finds per sequence.
        currentSize = currentSizes[i];
        foundMatches[i] = (int*) malloc(currentSize * sizeof(int));
        if (!foundMatches[i]) return 0;

        matchCounter[i] = 0;
        patternLength = patternLengths[i];
        currentPatterns = patterns[i];
        char commonStart[patternLength];
        if (patternCount[i] == 1) strcpy(commonStart, currentPatterns[0]);
        else {
            findCommonStart(commonStart, currentPatterns, patternLength,
                            i, patternCount);
        }

        if (threadsPerThread > 1) { // Parallel
            threadData.currIndex = i;
            threadData.currentSize = &currentSizes[i];
            threadData.currentPatterns = currentPatterns;
            threadData.commonStart = commonStart;
            size_t sequenceLength = strlen(sequences[i]);
            int portion = (int) sequenceLength / threadsPerThread;
            int remainder = (int) sequenceLength % threadsPerThread;
            for (int j = 0; j < threadsPerThread; j++) {
                if (j < remainder) {
                    threadData.portion = portion + 1;
                    threadData.start = j * threadData.portion;
                    threadData.end = threadData.start + threadData.portion;
                }
                else {
                    threadData.portion = portion;
                    threadData.start = portion*(j - remainder) + remainder*(portion + 1);
                    threadData.end = threadData.start + portion;
                }
                if (pthread_create(&threads[j], NULL, (void*) parallel_manageMatches, &threadData)) {
                    for (int k = 0; k < j; k++) {
                        pthread_join(threads[k], NULL);
                    }
                    free(threads);
                    pthread_barrier_destroy(&threadData.syncBarrier);
                    pthread_mutex_destroy(&threadData.resizeMutex);
                    pthread_mutex_destroy(&threadData.addMatchMutex);
                    return 0;
                }
                pthread_barrier_wait(&threadData.syncBarrier);
            }
            for (int j = 0; j < threadsPerThread; j++) {
                int result;
                pthread_join(threads[j], (void**) &result);
                if (!result) {
                    pthread_barrier_destroy(&threadData.syncBarrier);
                    pthread_mutex_destroy(&threadData.resizeMutex);
                    pthread_mutex_destroy(&threadData.addMatchMutex);
                    free(threads);
                    return 0;
                }
            }
        }
        else { // Sequential
            if (!findMatches(sequences[i], currentPatterns, commonStart, foundMatches,
                             &matchCounter[i], &currentSize, patternCount[i],
                             i, patternLength, 0, NULL, NULL,  0)) {
                return 0;
            }
        }
        foundMatches[i] = realloc(foundMatches[i], matchCounter[i] * sizeof(int));
        if (!foundMatches[i]) return 0;
    }
    if (threadsPerThread > 1) {
        pthread_barrier_destroy(&threadData.syncBarrier);
        pthread_mutex_destroy(&threadData.resizeMutex);
        pthread_mutex_destroy(&threadData.addMatchMutex);
        free(threads);
    }
    return 1;
}

int DNA_Search(FILE* file, int threadCount, int threadsPerThread) {
    if (threadCount < 1) return 0;
    int sequenceCount = 0;
    int maxSequenceLength = 0;

    if (!getSequenceData(file, &sequenceCount, &maxSequenceLength)) {
        return 0;
    }

    // Max 1 thread per sequence
    if (threadCount > sequenceCount) threadCount = sequenceCount;

    char*** patterns;
    size_t* patternLengths;
    if (!malloc_pattern_arrays(&patterns, &patternLengths, sequenceCount)) {
        return 0;
    }

    // Rough estimation that all patterns include 3 variations.
    size_t estimatedPatternCount = 3;
    if (!malloc_pattern_array(patterns,estimatedPatternCount, sequenceCount)) {
        return 0;
    }

    int patternCount[sequenceCount]; // Number of patterns for each sequence
    // Separate out pattern variations
    if (!distinguishPatterns(file, patterns, patternLengths, sequenceCount,
                             estimatedPatternCount, patternCount)) {
        return 0;
    }

    char** sequences = (char**) malloc(sequenceCount * sizeof(char*));
    if (!sequences) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        return 0;
    }

    // Fill "sequences" array with sequences found in the file.
    if (!parseSequences(file, sequences, sequenceCount, maxSequenceLength)) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        return 0;
    }

    pthread_t* threads = (pthread_t*) malloc(threadCount * sizeof(pthread_t));
    if (!threads) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        free_full_2DArray((void**) sequences, sequenceCount);
        return 0;
    }

    struct ThreadData data;
    if (pthread_barrier_init(&data.syncBarrier, NULL, 2)) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        free_full_2DArray((void**) sequences, sequenceCount);
        free(threads);
        return 0;
    }

    int** foundMatches = (int**) malloc(sequenceCount * sizeof(int*));
    if (!foundMatches) {
        free_pattern_arrays(patterns, patternLengths, sequenceCount,
                            sequenceCount, patternCount);
        free_full_2DArray((void**) sequences, sequenceCount);
        free(threads);
        pthread_barrier_destroy(&data.syncBarrier);
        return 0;
    }

    int matchCount[sequenceCount];
    struct SharedData sharedData = {
        .sequences = sequences,
        .patterns = patterns,
        .foundMatches = foundMatches,
        .patternLengths = patternLengths,
        .patternCount = patternCount,
        .matchCounter = matchCount,
    };
    data.sequenceCount = sequenceCount;
    data.sharedData = sharedData;
    data.extraThreads = threadsPerThread;

    int portion = sequenceCount / threadCount;
    int remainder = sequenceCount % threadCount;
    for (int i = 0; i < threadCount; i++) {
        if (i < remainder) {
            data.start = i * (portion + 1);
            data.end = data.start + portion + 1;
        }
        else {
            data.start = portion*(i - remainder) + remainder*(portion + 1);
            data.end = data.start + portion;
        }
        if (pthread_create(&threads[i], NULL, (void*) manageMatches, &data)) {
            printf("An error has occurred.\n");
            for (int j = 0; j < i; j++) {
                pthread_join(threads[j], NULL);
            }
            free_pattern_arrays(patterns, patternLengths, sequenceCount,
                                sequenceCount, patternCount);
            free_full_2DArray((void**) sequences, sequenceCount);
            free_full_2DArray((void**) foundMatches, sequenceCount);
            free(threads);
            pthread_barrier_destroy(&data.syncBarrier);
            return 0;
        }
        pthread_barrier_wait(&data.syncBarrier);
    }
    pthread_barrier_destroy(&data.syncBarrier);

    int result;
    for (int i = 0; i < threadCount; i++) {
        pthread_join(threads[i], (void**) &result);
        if (!result) {
            printf("An error has occurred.\n");
            free_pattern_arrays(patterns, patternLengths, sequenceCount,
                                sequenceCount, patternCount);
            free_full_2DArray((void **) sequences, sequenceCount);
            free(threads);
            return 0;
        }
    }
    free(threads);

    for (int i = 0; i < sequenceCount; i++) {
        printf("Occurrences in sequence %d: %d\n", i+1, matchCount[i]);
        for (int j = 0; j < matchCount[i]; j++) {
            printf("Occurrence %d at index: %d\n", j+1, foundMatches[i][j]);
        }
    }

    free_pattern_arrays(patterns, patternLengths, sequenceCount,
                        sequenceCount, patternCount);
    free_full_2DArray((void**) sequences, sequenceCount);
    free_full_2DArray((void**) foundMatches, sequenceCount);
    return 1;
}

int main() {
    FILE *file = getFile("sequences.txt");
    if (!file) {
        printf("File not found.\n");
        return -1;
    }
    if (!DNA_Search(file, 8, 3)) return -1;
    return 0;
}
