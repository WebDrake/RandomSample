// MyRandomSample
/**
Selects a random subsample out of $(D r), containing exactly $(D n)
elements. The order of elements is the same as in the original
range. The total length of $(D r) must be known. If $(D total) is
passed in, the total number of sample is considered to be $(D
total). Otherwise, $(D MyRandomSample) uses $(D r.length).

If the number of elements is not exactly $(D total), $(D
MyRandomSample) throws an exception. This is because $(D total) is
essential to computing the probability of selecting elements in the
range.

Example:
----
int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];
// Print 5 random elements picked off from a
foreach (e; myRandomSample(a, 5))
{
    writeln(e);
}
----
 */
struct MyRandomSample(R, Random = void)
    if(isUniformRNG!Random || is(Random == void))
{
    private size_t _available, _toSelect;
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
        _available = total;
        _toSelect = howMany;
        enforce(_toSelect <= _available);
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

    private void prime()
    {
        if (empty) return;
        assert(_available && _available >= _toSelect);
        for (;;)
        {
            static if(is(Random == void))
            {
                auto r = uniform(0, _available);
            }
            else
            {
                auto r = uniform(0, _available, gen);
            }

            if (r < _toSelect)
            {
                // chosen!
                return;
            }
            // not chosen, retry
            assert(!_input.empty);
            _input.popFront();
            ++_index;
            --_available;
            assert(_available > 0);
        }
    }
}

/// Ditto
auto myRandomSample(R)(R r, size_t n, size_t total)
if(isInputRange!R)
{
    return MyRandomSample!(R, void)(r, n, total);
}

/// Ditto
auto myRandomSample(R)(R r, size_t n) if (hasLength!R)
{
    return MyRandomSample!(R, void)(r, n, r.length);
}

/// Ditto
auto myRandomSample(R, Random)(R r, size_t n, size_t total, Random gen)
if(isInputRange!R && isUniformRNG!Random)
{
    auto ret = MyRandomSample!(R, Random)(r, n, total);
    ret.gen = gen;
    return ret;
}

/// Ditto
auto myRandomSample(R, Random)(R r, size_t n, Random gen)
if (isInputRange!R && hasLength!R && isUniformRNG!Random)
{
    auto ret = MyRandomSample!(R, Random)(r, n, r.length);
    ret.gen = gen;
    return ret;
}

unittest
{
    Random gen;
    int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];
    static assert(isForwardRange!(typeof(myRandomSample(a, 5))));
    static assert(isForwardRange!(typeof(myRandomSample(a, 5, gen))));

    //int[] a = [ 0, 1, 2 ];
    assert(myRandomSample(a, 5).length == 5);
    assert(myRandomSample(a, 5, 10).length == 5);
    assert(myRandomSample(a, 5, gen).length == 5);
    uint i;
    foreach (e; myRandomSample(randomCover(a, rndGen), 5))
    {
        ++i;
        //writeln(e);
    }
    assert(i == 5);
}

