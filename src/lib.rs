use error::Error;

pub mod error;

pub struct Stat
{
    dimension: usize,

    s1: usize,

    sv: Vec<f64>,
    svv: Vec<f64>,
    svvv: Vec<f64>,
    svvvv: Vec<f64>,

    sxy: Vec<Vec<f64>>,
    sxxy: Vec<Vec<f64>>,
    sxyy: Vec<Vec<f64>>,
    sxxyy: Vec<Vec<f64>>,
}

impl Stat
{
    pub fn new(dimension: usize) -> Self
    {
        Self
        {
            dimension,

            s1: 0,

            sv: vec![0.; dimension],
            svv: vec![0.; dimension],
            svvv: vec![0.; dimension],
            svvvv: vec![0.; dimension],

            sxy: vec![vec![0.; dimension]; dimension],
            sxxy: vec![vec![0.; dimension]; dimension],
            sxyy: vec![vec![0.; dimension]; dimension],
            sxxyy: vec![vec![0.; dimension]; dimension],
        }
    }

    pub fn fill(&mut self, input: &Vec<f64>) -> Result<(), Error>
    {
        if input.len() != self.dimension
        {
            return Err(Error::WrongInputSize)
        }

        self.s1 += 1;

        for (i, v_i) in input.iter().enumerate()
        {
            self.sv[i] += v_i;
            self.svv[i] += v_i*v_i;
            self.svvv[i] += v_i*v_i*v_i;
            self.svvvv[i] += v_i*v_i*v_i*v_i;
        
            for (j, v_j) in input.iter().enumerate()
            {
				self.sxy[i][j] += v_i * v_j;
				self.sxxy[i][j] += v_i*v_i * v_j;
				self.sxyy[i][j] += v_i * v_j*v_j;
				self.sxxyy[i][j] += v_i*v_i * v_j*v_j;
            }
        }

        Ok(())
    }

    /// Returns the number of entries.
    pub fn get_entries(&self) -> usize
    {
        self.s1
    }

    fn assert_index_in_range(&self, idx: usize) -> Result<(), Error>
    {
        if idx >= self.dimension
        {
            return Err(Error::IndexOutOfRange)
        }

        Ok(())
    }

    fn assert_size_at_least(&self, min_size: usize) -> Result<(), Error>
    {
        if self.s1 < min_size
        {
            return Err(Error::SampleTooSmall)
        }

        Ok(())
    }
    
    // -------------------- 1D getterns --------------------------------------------------

    /// Returns the mean or None if empty.
    pub fn get_mean(&self, idx: usize) -> Result<f64, Error>
    {
        self.assert_index_in_range(idx)?;
        self.assert_size_at_least(1)?;
        Ok(self.sv[idx] / (self.s1 as f64))
    }

    /// Returns the mean uncertainty.
    pub fn get_mean_unc(&self, idx: usize) -> Result<f64, Error>
    {
        Ok(self.get_std_dev(idx)? / (self.s1 as f64).sqrt())
    }

    /// Returns the standard deviation.
    pub fn get_std_dev(&self, idx: usize) -> Result<f64, Error>
    {
        self.assert_index_in_range(idx)?;
        self.assert_size_at_least(2)?;

        let s1 = self.s1 as f64;
        let v = (self.svv[idx] - self.sv[idx]*self.sv[idx] / s1) / (s1 - 1.);
        Ok(if v > 0. { v.sqrt() } else { 0. })
    }

    /// Returns the standard deviation uncertainty.
    pub fn get_std_dev_unc(&self, idx: usize) -> Result<f64, Error>
    {
        self.assert_index_in_range(idx)?;
        self.assert_size_at_least(4)?;

        let mu = self.get_mean(idx)?;
        let s = self.get_std_dev(idx)?;
        let v = s * s;
    
        let s1 = self.s1 as f64;
        let sum = self.svvvv[idx] - 4.*mu*self.svvv[idx] + 6.*mu*mu*self.svv[idx] - 4.*mu*mu*mu*self.sv[idx] + mu*mu*mu*mu*s1;
        let e4 = sum / (s1 - 1.);
    
        let v_var = (e4 - (s1 - 3.) / (s1 - 1.)*v*v) / s1;
        let s_var = v_var / 4. / v;
        let s_s = if s_var > 0. { s_var.sqrt() } else { 0. };
        Ok(s_s)
    }

    /// Returns an approximation of the standard deviation uncertainty valid for Gaussian populations.
    pub fn get_std_dev_unc_gauss(&self, idx: usize) -> Result<f64, Error>
    {
        let s = self.get_std_dev(idx)?;
	    Ok(s / (2. * (self.s1 as f64)).sqrt())
    }
    
    // -------------------- 2D getterns --------------------------------------------------

    pub fn get_covariance(&self, i: usize, j: usize) -> Result<f64, Error>
    {
        self.assert_index_in_range(i)?;
        self.assert_index_in_range(j)?;
        self.assert_size_at_least(2)?;

        let s1 = self.s1 as f64;
        Ok((self.sxy[i][j] - self.sv[i]*self.sv[j] / s1) / (s1 - 1.))
    }

	pub fn get_correlation(&self, i: usize, j: usize) -> Result<f64, Error>
    {
        let num = self.get_covariance(i, j)?;
        let den = self.get_std_dev(i)? * self.get_std_dev(j)?;
        if den <= 0.
        {
            return Err(Error::ZeroDenominator);
        }

        Ok(num / den)
    }

	pub fn get_covariance_unc(&self, _i: usize, _j: usize) -> Result<f64, Error>
    {
        todo!();

        // TODO:
        /*
        double mx = GetMean(i);
        double my = GetMean(j);
        double sx = GetStdDev(i);
        double sy = GetStdDev(j);
        double C = GetCovariance(i, j);

        double sum =
            Sxxyy[i][j] 
            -2.*Sxyy[i][j]*mx - 2.*Sxxy[i][j]*my
            + 4.*Sxy[i][j]*mx*my
            + Svv[i]*my*my + Svv[j]*mx*mx
            - 2.*Sv[i]*mx*my*my - 2.*Sv[j]*mx*mx*my
            + mx*mx*my*my;
        double D = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

        double C_var = (S1 > 2.) ? (D + sx*sx*sy*sy/(S1 - 1.) - (S1-2.)/(S1-1.)*C*C) / S1 : 0.;
        double C_s = (C_var > 0.) ? sqrt(C_var) : 0.;

        return C_s;
        */
    }

	pub fn get_correlation_unc(&self, _i: usize, _j: usize) -> Result<f64, Error>
    {
        todo!();

        // TODO:
        // WARNING: the calculation below assumes no correlation between C, si_i and si_j, which
        // might not be correct - in that case it gives an upper bound for the uncertainty

        /*
        double C = GetCovariance(i, j), C_unc = GetCovarianceUnc(i, j);
        double si_i = GetStdDev(i), si_i_unc = GetStdDevUnc(i);
        double si_j = GetStdDev(j), si_j_unc = GetStdDevUnc(j);
        double rho = C / (si_i * si_j);
        double sum =
            (C != 0. && si_i != 0. && si_j != 0.) ? pow(C_unc / C, 2.) + pow(si_i_unc / si_i, 2.) + pow(si_j_unc / si_j, 2.) : 0.;
        double rho_unc = fabs(rho) * sqrt(sum);
        return rho_unc;
        */
    }

    // -------------------- print methods --------------------------------------------------

	pub fn print_stat(&self)
    {
        println!("entries: {:.3E}", self.s1);
    }

	pub fn print_mean_and_std_dev(&self) -> Result<(), Error>
    {
        for i in 0..self.dimension
        {
            println!("{}: mean {:+.3E} +- {:.3E}, std. dev. = {:.3E} +- {:.3E}", i,
                self.get_mean(i)?, self.get_mean_unc(i)?,
                self.get_std_dev(i)?, self.get_std_dev_unc(i)?);
        }

        Ok(())
    }
	
	pub fn print_covariance(&self) -> Result<(), Error>
    {
        print!("      ");
        for i in 0..self.dimension
        {
            print!("   {:6}", i);
        }
        println!("");

        for i in 0..self.dimension
        {
            print!("{:6}", i);

            for j in 0..self.dimension
            {
                print!("   {:+.3}", self.get_covariance(i, j)?);
            }

            println!();
        }

        Ok(())
    }

	pub fn print_correlation(&self) -> Result<(), Error>
    {
        print!("      ");
        for i in 0..self.dimension
        {
            print!("   {:6}", i);
        }
        println!("");

        for i in 0..self.dimension
        {
            print!("{:6}", i);

            for j in 0..self.dimension
            {
                print!("   {:+.3}", self.get_correlation(i, j)?);
            }

            println!();
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    fn assert_error_eq<T>(input: Result<T, Error>, expected_error: Error)
    {
        match input
        {
            Ok(_) => panic!("result is not error"),
            Err(e) => {
                if e != expected_error
                {
                    panic!("the error is not the expected one")
                }
            },
        }
    }

    fn assert_eq_with_tolerance(left: f64, right: f64, tolerance: f64)
    {
        if (left - right).abs() > tolerance
        {
            panic!("{} is incompatible with {} given tolerance {}", left, right, tolerance)
        }
    }

    #[test]
    fn test_out_of_range()
    {
        let s = Stat::new(1);

        assert_error_eq(s.get_mean(17), Error::IndexOutOfRange);
        assert_error_eq(s.get_mean_unc(17), Error::IndexOutOfRange);
        assert_error_eq(s.get_std_dev(17), Error::IndexOutOfRange);
        assert_error_eq(s.get_std_dev_unc(17), Error::IndexOutOfRange);
        assert_error_eq(s.get_std_dev_unc_gauss(17), Error::IndexOutOfRange);
    }

    #[test]
    fn test_empty()
    {
        let s = Stat::new(1);

        assert_error_eq(s.get_mean(0), Error::SampleTooSmall);
        assert_error_eq(s.get_mean_unc(0), Error::SampleTooSmall);
        assert_error_eq(s.get_std_dev(0), Error::SampleTooSmall);
        assert_error_eq(s.get_std_dev_unc(0), Error::SampleTooSmall);
        assert_error_eq(s.get_std_dev_unc_gauss(0), Error::SampleTooSmall);
    }

    fn build_std_stat() -> Stat
    {
        let mut s = Stat::new(2);

        s.fill(& vec![1., 4.]).unwrap();
        s.fill(& vec![3., 0.]).unwrap();
        s.fill(& vec![2., -1.]).unwrap();
        s.fill(& vec![0., 2.]).unwrap();

        s
    }

    #[test]
    fn test_std()
    {
        let s = build_std_stat();

        s.print_stat();
        s.print_mean_and_std_dev().unwrap();
        s.print_covariance().unwrap();
        s.print_correlation().unwrap();

        assert_eq_with_tolerance(s.get_mean(0).unwrap(), 1.5, 1E-6);
        assert_eq_with_tolerance(s.get_mean(1).unwrap(), 1.25, 1E-6);
        
        assert_eq_with_tolerance(s.get_mean_unc(0).unwrap(), 0.64549722, 1E-6);

        assert_eq_with_tolerance(s.get_std_dev(0).unwrap(), 1.2909944, 1E-6);

        assert_eq_with_tolerance(s.get_std_dev_unc_gauss(0).unwrap(), 0.45643546, 1E-6);

        // TODO: test correlations and covariances
    }
}