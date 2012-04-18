import std.algorithm, std.c.time, std.conv, std.exception,
       std.math, std.numeric, std.range, std.traits,
       core.thread, core.time;

import std.random, std.stdio;

// RandomSampleVitter
/**
Selects a random subsample out of $(D r), containing exactly $(D n)
elements. The order of elements is the same as in the original
range. The total length of $(D r) must be known. If $(D total) is
passed in, the total number of sample is considered to be $(D
total). Otherwise, $(D RandomSampleVitter) uses $(D r.length).

If the number of elements is not exactly $(D total), $(D
RandomSampleVitter) throws an exception. This is because $(D total) is
essential to computing the probability of selecting elements in the
range.

$(D RandomSampleVitter) implements Jeffrey Scott Vitter's Algorithm D,
which selects a sample of size $(D n) in O(n) steps and requiring O(n)
random variates, regardless of the size of the data being sampled.

Example:
----
int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];
// Print 5 random elements picked off from a
foreach (e; randomSampleVitter(a, 5))
{
    writeln(e);
}
----
 */
struct RandomSampleVitter(R, Random = void)
    if(isUniformRNG!Random || is(Random == void))
{
    private size_t _available, _howMany, _toSelect, _total;
    private immutable ushort _alphaInverse = 13; // Vitter's recommended value.
    private bool _algorithmA;
    private double _Vprime;
    private R _input;
    private size_t _index;

    // If we're using the default thread-local random number generator then
    // we shouldn't store a copy of it here.  Random == void is a sentinel
    // for this.  If we're using a user-specified generator then we have no
    // choice but to store a copy.
    static if(!is(Random == void))
    {
        Random gen;
    }

/**
Constructor.
*/
    static if (hasLength!R)
        this(R input, size_t howMany)
        {
            this(input, howMany, input.length);
        }

    this(R input, size_t howMany, size_t total)
    {
        _input = input;
        _available = _total = total;
        _toSelect = _howMany = howMany;
        enforce(_toSelect <= _available);

        // We can save ourselves a random variate by checking right
        // at the beginning if we should use Algorithm A.
        if((_alphaInverse * _toSelect) > _available)
        {
            _algorithmA = true;
        }
        else
        {
            newVprime(_toSelect);
            _algorithmA = false;
        }
        // we should skip some elements initially so we don't always
        // start with the first
        prime();
    }

/**
   Range primitives.
*/
    @property bool empty() const
    {
        return _toSelect == 0;
    }

    @property auto ref front()
    {
        assert(!empty);
        return _input.front;
    }

/// Ditto
    void popFront()
    {
        _input.popFront();
        --_available;
        --_toSelect;
        ++_index;
        prime();
    }

/// Ditto
    @property typeof(this) save()
    {
        auto ret = this;
        ret._input = _input.save;
        return ret;
    }

/// Ditto
    @property size_t length()
    {
        return _toSelect;
    }

/**
Returns the index of the visited record.
 */
    size_t index()
    {
        return _index;
    }

/**
Vitter's Algorithm A, used when the ratio of needed sample values
to remaining data values is sufficiently large.
*/
    private size_t skipA()
    {
        size_t S;
        double V, quot, top;

        if(_toSelect==1)
        {
            static if(is(Random==void))
            {
                S = uniform(0, _available);
            }
            else
            {
                S = uniform(0, _available, gen);
            }
        }
        else {
            S = 0;
            top = _available - _toSelect;
            quot = top / _available;

            static if(is(Random==void))
            {
                V = uniform!("()")(0.0, 1.0);
            }
            else
            {
                V = uniform!("()")(0.0, 1.0, gen);
            }

            while (quot > V) {
                ++S;
                quot *= (top - S) / (_available - S);
            }
        }

        return S;
    }

/**
Randomly reset the value of _Vprime.
*/
    private void newVprime(size_t remaining)
    {
        static if(is(Random == void))
        {
            double r = uniform!("()")(0.0, 1.0);
        }
        else
        {
            double r = uniform!("()")(0.0, 1.0, gen);
        }

        _Vprime = r ^^ (1.0 / remaining);
    }

/**
Vitter's Algorithm D.  For an extensive description of the algorithm
and its rationale, see:

  * Vitter, J.S. (1984), "Faster methods for random sampling",
    Commun. ACM 27(7): 703--718

  * Vitter, J.S. (1987) "An efficient algorithm for sequential random
    sampling", ACM Trans. Math. Softw. 13(1): 58-67.

Variable names are chosen to match those in Vitter's paper.
*/
    private size_t skip()
    {
        // If the number of points still to select is greater than a
        // certain proportion of the remaining data points, we carry
        // out the sampling with Algorithm A.
        if(_algorithmA)
        {
            return skipA;
        }
        else if((_alphaInverse * _toSelect) > _available)
        {
            _algorithmA = true;
            return skipA;
        }
        // Otherwise, we use the standard Algorithm D mechanism.
        else if ( _toSelect > 1 )
        {
            size_t S;
            size_t qu1 = 1 + _available - _toSelect;
            double X, y1;

            while(1)
            {
                // Step D2: set values of X and U.
                for(X = _available * (1-_Vprime), S = cast(size_t) trunc(X);
                    S >= qu1;
                    X = _available * (1-_Vprime), S = cast(size_t) trunc(X))
                {
                    newVprime(_toSelect);
                }

                static if(is(Random == void))
                {
                    double U = uniform!("()")(0.0, 1.0);
                }
                else
                {
                    double U = uniform!("()")(0.0, 1.0, gen);
                }

                y1 = (U * (cast(double) _available) / qu1) ^^ (1.0/(_toSelect - 1));

                _Vprime = y1 * ((-X/_available)+1.0) * ( qu1/( (cast(double) qu1) - S ) );

                // Step D3: if _Vprime <= 1.0 our work is done and we return S.
                // Otherwise ...
                if(_Vprime > 1.0)
                {
                    size_t top = _available - 1, limit;
                    double y2 = 1.0, bottom;

                    if(_toSelect > (S+1) )
                    {
                        bottom = _available - _toSelect;
                        limit = _available - S;
                    }
                    else
                    {
                        bottom = _available - (S+1);
                        limit = qu1;
                    }

                    foreach(size_t t; limit.._available)
                    {
                        y2 *= top/bottom;
                        top--;
                        bottom--;
                    }

                    // Step D4: decide whether or not to accept the current value of S.
                    if( (_available/(_available-X)) < (y1 * (y2 ^^ (1.0/(_toSelect-1)))) )
                    {
                        // If it's not acceptable, we generate a new value of _Vprime
                        // and go back to the start of the for(;;) loop.
                        newVprime(_toSelect);
                    }
                    else
                    {
                        // If it's acceptable we generate a new value of _Vprime
                        // based on the remaining number of sample points needed,
                        // and return S.
                        newVprime(_toSelect-1);
                        return S;
                    }
                }
                else
                {
                    // Return if condition D3 satisfied.
                    return S;
                }
            }
        }
        else
        {
            // If only one sample point remains to be taken ...
            return cast(size_t) trunc(_available * _Vprime);
        }
    }

    private void prime()
    {
        if (empty) return;
        assert(_available && _available >= _toSelect);
        immutable size_t S = skip;
        _input.popFrontN(S);
        _index += S;
        _available -= S;
        assert(_available > 0);
        return;
    }
}

/// Ditto
auto randomSampleVitter(R)(R r, size_t n, size_t total)
if(isInputRange!R)
{
    return RandomSampleVitter!(R, void)(r, n, total);
}

/// Ditto
auto randomSampleVitter(R)(R r, size_t n) if (hasLength!R)
{
    return RandomSampleVitter!(R, void)(r, n, r.length);
}

/// Ditto
auto randomSampleVitter(R, Random)(R r, size_t n, size_t total, Random gen)
if(isInputRange!R && isUniformRNG!Random)
{
    auto ret = RandomSampleVitter!(R, Random)(r, n, total);
    ret.gen = gen;
    return ret;
}

/// Ditto
auto randomSampleVitter(R, Random)(R r, size_t n, Random gen)
if (isInputRange!R && hasLength!R && isUniformRNG!Random)
{
    auto ret = RandomSampleVitter!(R, Random)(r, n, r.length);
    ret.gen = gen;
    return ret;
}

unittest
{
    Random gen;
    int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];
    static assert(isForwardRange!(typeof(randomSampleVitter(a, 5))));
    static assert(isForwardRange!(typeof(randomSampleVitter(a, 5, gen))));

    //int[] a = [ 0, 1, 2 ];
    assert(randomSampleVitter(a, 5).length == 5);
    assert(randomSampleVitter(a, 5, 10).length == 5);
    assert(randomSampleVitter(a, 5, gen).length == 5);
    uint i;
    foreach (e; randomSampleVitter(randomCover(a, rndGen), 5))
    {
        ++i;
        //writeln(e);
    }
    assert(i == 5);
}

void samplingTestAggregate(size_t total, size_t n, size_t repeats=1, bool verbose=true)
{
    double[] recordCountS, recordCountD;
    clock_t start_time, end_time;

    void displayResults(double[] recordCount)
    {
        foreach(size_t i, double c; recordCount)
            writeln("\trecord ", i, " was picked ", c, " times.");
    }

    writeln("Picking ", n, " from ", total, ", ", repeats, " times.");
    writeln;

    writeln("Algorithm S:");
    recordCountS.length = total;
    recordCountS[] = 0.0;
    start_time = clock();
    foreach(size_t i; 0..repeats)
    {
        auto sampleS = randomSample(iota(0, total), n);
        foreach(size_t s; sampleS)
            recordCountS[s]++;
    }
    end_time = clock();
    if(verbose)
        displayResults(recordCountS);
    writeln("\t\tSampling completed in ",(cast(double) (end_time - start_time))/CLOCKS_PER_SEC, " seconds.");
    writeln;

    writeln("Algorithm D:");
    recordCountD.length = total;
    recordCountD[] = 0.0;
    start_time = clock();
    foreach(size_t i; 0..repeats)
    {
        auto sampleD = randomSampleVitter(iota(0, total), n);
        foreach(size_t s; sampleD)
            recordCountD[s]++;
    }
    end_time = clock();
    if(verbose)
        displayResults(recordCountD);
    writeln("\t\tSampling completed in ",(cast(double) (end_time - start_time))/CLOCKS_PER_SEC, " seconds.");
    writeln;
}

void main(string[] args)
{
    auto s = randomSampleVitter(iota(0,100),5);

    foreach(uint i; s)
        writeln(i);

    writeln();

    foreach(uint i; s)
        writeln(i);

    samplingTestAggregate(100, 5, 1_000_000);

    samplingTestAggregate(1000, 5, 1_000_000, false);

    samplingTestAggregate(100_000_000, 1_000, 1, false);
}

