# breakage-points
Breakage Point Homework Assignment

## Objectives

1. Parse the Tab separated file. Each line specifies a "translocation".

2. Figure out the two most likely FragileSites that yielded a given
Translocation and the coordinates of these Fragile sites. It is not required to do this for the entire dataset.

3. Serialize the resulting set of Fragile Sites to something easy to parse.
YAML (via PyYaml) and CSV (document your field names) via the built in
lib are both good, basic choices.

## General Advice

There are lots of ways to do this, and I want to see **your**
approach, but here is a hint from the DB authors:
'''
More advanced users might prefer to run customized BLAST searches and
inspect the output in order to identify the precise breakpoints, which
is the process we followed in the original BMC Genomics paper. This
has the advantage that it will uncover microhomologies present at the
fusion boundary (there are quite a few of those).
'''
I want to see your choice and **why**. For example "Used BLAST API
because I already have a good client for it" or "Wrote out FASTQ file
and then invoked bwa-mem as a subprocess because it is fast"

__Bonus Points - de-duplicate the set of resulting FragileSites based
on their coordinates. That is: if two FragileSites have the same
genomic coordinates, then they must be the same fragile site__

Equally, you should demonstrate **how you know** you've found
reasonable genomic coordinates. This means you should write an
automated test. For example you could assert that the coordinates
occur in the stated gene for a sample of results.

It is **not** important to do this for the entire data set. I don't
actually care about the actual output or even speed. I care about your
design judgement, if you write clean, PEP08 compliant Python, how you
tested it to show it both runs and produces reasonably correct output,
and how you document and comment it so that other team members can
understand.

If you're not up to date on professional Python coding style check out
some popular libraries, such as SQLAlchemy, or "Code Like a
Pythonista" 

See http://python.net/~goodger/projects/pycon/2007/idiomatic/handout.html

Code quality
([readability](http://python.net/~goodger/projects/pycon/2007/idiomatic/handout.html#coding-style-readability-counts),
consistency) is not emphasised much in academia, but a lack of it drives
professional developers **insane** (in a  way that is probably unhealthy).
Most developers aren't biologists and will have no idea if your implementation
is biologically accurate (that is your concern!), but they will be able to
assess its logic and consistency. Issues with code style and readability will
present a roadblock to that.

ALWAYS make sure you use the `if __name__ == '__main__'` idiom in your
scripts. I can't tell you how many candidates forget this, causing
their code to execute on import.

There are many tools such as PyLint that will check this for you
automatically. Most IDEs have this built in. PyCharm has a free
community edition and a free or steeply discounted Educational edition
that works well.

## Data Model
__We shall define a Translocation object as a type of GenomeFeature
composed of 2 FragileSite objects (which are a subclass of
GenomeFeature)  such that `fragile_site_a + fragile_site_b` -->
`translocation_a_b`__

A GenomeFeature is an object that has at least
 1) a name and
 2) various other attributes and
 3) coordinates (in a genome)

Coordinates are a tuple of `(chromosome, absolute_start, absolute_end,
strand)`. For example `(1, 123, 456, -1)` would be a Feature that occurs on
Chromosome `1` from `123` to `456` on the negative strand, such that it's
sequence is `reverse_compliment(chromosome_one.sequence[123:456])`, that is you
can use Python string slicing to get the positive strand sequence, but need to
reverse compliment that sequence if the feature is on the negative strand to
get its sequence 5' to 3'

## Scaffold
A basic scaffold for a Python project has been provided in this repository. It
is recommended that a
[virtualenv](http://python.net/~goodger/projects/pycon/2007/idiomatic/handout.html#coding-style-readability-counts)
is used for your development environment. The project can be installed in
"development" mode by running the following command (whilst the virtualenv is
activated) in the root of the repository:

``` bash
$ python setup.py develop
```

This will install the module `brkpoints` and any dependencies listed in the
`setup.py` file, after this the test suite can be run with the following
command:

``` bash
$ nosetests
```

To see more verbose output and the output from `print` statements, invoke with the following flags:

``` bash
$ nosetests -v --nocapture
```
