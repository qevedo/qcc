//! The Hungarian algorithm: the cheapest way to give each row a different
//! column of a cost matrix.

/// Costs at or above this are "forbidden" (a full core, for example).
pub const FORBIDDEN: f64 = 1e12;

/// For a rectangular `cost` matrix (rows × columns), assign as many rows as
/// possible to distinct columns at the minimum total cost: `result[row]` is
/// the column of `row`, or `None` when it got none (more rows than columns,
/// or only forbidden columns left).
///
/// Kuhn–Munkres with potentials, O(n² m) for n ≤ m.
pub fn assign(cost: &[Vec<f64>]) -> Vec<Option<usize>> {
    let rows = cost.len();
    if rows == 0 {
        return vec![];
    }
    let cols = cost[0].len();
    if cols == 0 {
        return vec![None; rows];
    }
    // The algorithm needs rows ≤ columns: solve the transpose otherwise.
    let mut result = if rows <= cols {
        solve(cost)
    } else {
        let transposed: Vec<Vec<f64>> = (0..cols)
            .map(|c| (0..rows).map(|r| cost[r][c]).collect())
            .collect();
        let mut by_row = vec![None; rows];
        for (c, r) in solve(&transposed).into_iter().enumerate() {
            if let Some(r) = r {
                by_row[r] = Some(c);
            }
        }
        by_row
    };
    for (r, c) in result.iter_mut().enumerate() {
        if c.is_some_and(|c| cost[r][c] >= FORBIDDEN) {
            *c = None;
        }
    }
    result
}

/// The square-or-wide case (rows ≤ columns), 1-indexed internally.
fn solve(a: &[Vec<f64>]) -> Vec<Option<usize>> {
    let (n, m) = (a.len(), a[0].len());
    let cost = |i: usize, j: usize| a[i - 1][j - 1].min(FORBIDDEN);
    let mut u = vec![0.0; n + 1];
    let mut v = vec![0.0; m + 1];
    let mut p = vec![0usize; m + 1]; // p[j]: the row matched to column j
    let mut way = vec![0usize; m + 1];
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut minv = vec![f64::INFINITY; m + 1];
        let mut used = vec![false; m + 1];
        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = f64::INFINITY;
            let mut j1 = 0;
            for j in 1..=m {
                if !used[j] {
                    let cur = cost(i0, j) - u[i0] - v[j];
                    if cur < minv[j] {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if minv[j] < delta {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for j in 0..=m {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut result = vec![None; n];
    for j in 1..=m {
        if p[j] != 0 {
            result[p[j] - 1] = Some(j - 1);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn total(cost: &[Vec<f64>], a: &[Option<usize>]) -> f64 {
        a.iter()
            .enumerate()
            .filter_map(|(r, c)| c.map(|c| cost[r][c]))
            .sum()
    }

    /// The cheapest assignment by trying every permutation (small cases).
    fn brute_force(cost: &[Vec<f64>]) -> f64 {
        fn go(cost: &[Vec<f64>], row: usize, used: &mut Vec<bool>) -> f64 {
            if row == cost.len() {
                return 0.0;
            }
            let mut best = f64::INFINITY;
            for c in 0..used.len() {
                if !used[c] {
                    used[c] = true;
                    best = best.min(cost[row][c] + go(cost, row + 1, used));
                    used[c] = false;
                }
            }
            best
        }
        go(cost, 0, &mut vec![false; cost[0].len()])
    }

    #[test]
    fn matches_brute_force() {
        let mut rng = crate::rng::Rng::new(3);
        for _ in 0..200 {
            let (rows, cols) = (1 + rng.below(5), 5 + rng.below(2));
            let cost: Vec<Vec<f64>> = (0..rows)
                .map(|_| (0..cols).map(|_| rng.range(0.0, 2.0)).collect())
                .collect();
            let a = assign(&cost);
            assert!((total(&cost, &a) - brute_force(&cost)).abs() < 1e-9);
            let mut columns: Vec<usize> = a.iter().flatten().copied().collect();
            columns.sort();
            columns.dedup();
            assert_eq!(columns.len(), rows);
        }
    }

    #[test]
    fn more_rows_than_columns_and_forbidden_cells() {
        // Three operations, two cores, the second core full.
        let cost = vec![
            vec![1.0, FORBIDDEN],
            vec![0.5, FORBIDDEN],
            vec![2.0, FORBIDDEN],
        ];
        let a = assign(&cost);
        assert_eq!(a, vec![None, Some(0), None]);
    }
}
