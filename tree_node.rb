require 'matrix'

# Nested Interval Tree Encoding with Matrices
#
# Encodes a tree similar to a Materialized Path encoding:
#
#    food       = TreeNode[1]         1
#    meat       = TreeNode[1,1]       ├ 1.1
#    chicken    = TreeNode[1,1,1]     │  ├ 1.1.1
#    beef       = TreeNode[1,1,2]     │  ├ 1.1.2
#    steak      = TreeNode[1,1,2,1]   │  │  └ 1.1.2.1
#    pork       = TreeNode[1,1,3]     │  └ 1.1.3
#    dairy      = TreeNode[1,2]       ├ 1.2
#    milk       = TreeNode[1,2,1]     │  └ 1.2.1
#    produce    = TreeNode[1,3]       └ 1.3
#    fruit      = TreeNode[1,3,1]        ├ 1.3.1
#    vegetables = TreeNode[1,3,2]        └ 1.3.2
#    vehicle    = TreeNode[2]         2
#    automobile = TreeNode[2,1]       ├ 2.1
#    airplane   = TreeNode[2,2]       └ 2.2
#
# A node's interval encompasses all its descendants as in Nested Interval encoding:
#
#    f = food.interval
#    m = milk.interval
#    f.first <= m.first && m.first < m.last && m.last < m.last # => true
#    food.descendant?(milk) # more convenient
#
# Everything is encoded in a node's uniquely identifying 'id', a pair of integers:
#
#    milk.id # => [2, 1]
#    milk.id # => [8, 5]
#
# For RDBMS storage, ids and intervals need to be stored and indexed
#
# Acknowledgements: Vadim Tropashko
# https://vadimtropashko.files.wordpress.com/2011/07/ch5.pdf
class MatrixNestedInterval::TreeNode
  class << self
    def [](*path)
      all_positive_integers = path.all? { |n| n.integer? && n > 0 }
      raise ArgumentError, "path must consist only of positive integers" unless all_positive_integers
      matrix = path.inject(Matrix.identity(2)) { |mat, n| mat * atomic_matrix(n) }
      new(matrix)
    end

    def from_id(x, y)
      s, t = extended_euclidean(x, y)
      new(Matrix[[x, -t], [y, s]])
    end

    private

    def atomic_matrix(n)
      Matrix[[n.to_i + 1, -1], [1,  0]]
    end

    # Extended Euclidean algorithm
    #
    # given integers a and m, iteratively computes:
    #  * coefficients s and t that satisfies Bézout's Identity: gcd(a, m) == s*a + t*m
    #  * greatest common denominator: gcd(a, m)
    #  * modular multiplicative inverse of a with respect to the modulus m, if it exists
    def extended_euclidean(a, m)
      r0, r1 = a, m
      s, t = 1, 0
      until r1.zero?
        q, r2 = r0.divmod(r1) # Euclidean algorithm
        r0, r1 = r1, r2
        s, t = t, s - q * t # extended Euclidean algorithm
      end
      gcd = r0
      t = (gcd - s * a) / m # from Bézout's Identity
      inv = if gcd.abs == 1
              s < 0 ? s + m : s # modular multiplicative inverse
            end
      [s, t, gcd, inv]
    end
  end

  attr_reader :matrix

  def initialize(matrix)
    @matrix = matrix
    raise ArgumentError, "invalid #{self.class.name} id: #{id}" unless valid?
  end

  def id
    [a, c]
  end

  def ==(other)
    other.is_a?(self.class) && id == other.id
  end

  def n
    -1 - a.div(b) # NOTE: intentional truncation
  end

  def root
    return self if root?
    root_n = -1 - a.div(-c) # NOTE: intentional truncation
    self.class[root_n]
  end

  def root?
    d.zero?
  end

  def parent
    return nil if root?
    self.class.new(Matrix[[-b, a % b],
                          [-d, c % d]])
  end

  def child?(child)
    a == -child.b && c == -child.d
  end

  def interval
    left = Rational(a + b, c + d)
    right = Rational(a, c)
    left...right
  end

  def descendant?(other)
    self != other && interval.include?(other.interval)
  end

  def ancestor?(other)
    other != self && other.interval.include?(interval)
  end

  def path
    lineage.map(&:n)
  end

  def lineage
    Enumerator.new do |y|
      ancestor = self.root
      while ancestor
        y.yield(ancestor)
        cutting = cutting_from(ancestor)
        ancestor = cutting && cutting.root.grafted_onto(ancestor)
      end
    end
  end

  def ancestors
    Enumerator.new do |y|
      node = self
      y.yield(node) while node = node.parent
    end
  end

  def grafted_onto(other)
    self.class.new(other.matrix * @matrix)
  end

  def cutting_from(ancestor)
    return nil unless ancestor?(ancestor)
    inverse_ancestor = ancestor.matrix.inverse.map(&:to_i)
    self.class.new(inverse_ancestor * @matrix)
  end

  private

  def a; @matrix[0, 0]; end
  def b; @matrix[0, 1]; end
  def c; @matrix[1, 0]; end
  def d; @matrix[1, 1]; end

  def valid?
    (
      @matrix.row_size == 2 && @matrix.column_size == 2 &&
      (1...a).cover?(c) && (1...a).cover?(-b) &&
      a.gcd(c) == 1 && b.gcd(d) == 1 &&
      @matrix.determinant == 1
    )
  end
end
