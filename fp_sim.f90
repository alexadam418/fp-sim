
!--------------------------------------------------------------
!This module is to be used for simulating a fabry-perot spectrometer
!
!                 THE JANK IS DANK
!         
!               written by Alex Adam
!---------------------------------------------------------------

module fp_sim
      implicit none
      private
      double precision, public, parameter :: c = 299792458.
      double precision, public, parameter :: pi = 3.14159265359, &
                                              e = 2.71828, &
                                             kb = 1.380649D-23 , &
                                              h = 6.62607015D-34
      public  :: time_round_trip, &
                 photon_decay_time, &
                 cw_decay_rate, &
                 phase_shift, &
                 mode_index, &
                 free_spectral_range, &
                 wavenumber, &
                 lorentzian_FWHM, &
                 spectral_line_shape, &
                 lorentzian_lines, &
                 resonant_freq, &
                 create_spectrum, &
                 planck_spectrum, &
                 add_noise, &
                 scan_cavity

                              
     contains
     function time_round_trip(length) result(time_result)
          double precision, intent(in) :: length
          double precision :: time_result
          time_result = (2*length)/c
     end function time_round_trip

     function photon_decay_time(R_1, R_2, t_rt) result(time_const)
          double precision, intent(in) :: R_1, R_2, t_rt     
          double precision :: time_const
          time_const = t_rt/((-1*log(R_1)) + (-1*log(R_2)))
     end function photon_decay_time 
     
     function cw_decay_rate(tau_c, no_photons) result(rate)
          double precision, intent(in) :: tau_c, no_photons
          double precision :: rate
          rate = (1/tau_c)*no_photons
     end function cw_decay_rate

     function phase_shift(freq, t_rt) result(phase)
          double precision, intent(in) :: freq, t_rt
          double precision :: phase
          phase = pi*freq*t_rt
     end function phase_shift

     function mode_index(FSR, res_freq) result(q)
          double precision, intent(in) :: FSR, res_freq
          double precision :: q
          q = res_freq/FSR
     end function mode_index
      
     function free_spectral_range(t_rt) result(fs_range)
          double precision, intent(in) :: t_rt
          double precision :: fs_range
          fs_range = 1/t_rt
     end function
     
     function wavenumber(q, FSR) result(k)
          double precision, intent(in) :: q, FSR
          double precision :: k
          k = (2*pi*q*FSR)/c
     end function wavenumber

     function lorentzian_FWHM(tau_c) result(FWHM)
          double precision, intent(in) :: tau_c 
          double precision :: FWHM
          FWHM = 1/(2*pi*tau_c)
     end function lorentzian_FWHM

     function spectral_line_shape(res_frq, freq, FWHM) result(gamma)
          double precision, intent(in) :: res_frq, freq, FWHM
          double precision :: gamma
          gamma = (1/pi)*(FWHM/2)/((FWHM/2)**2 + (freq - res_frq)**2)
     end function spectral_line_shape

     function lorentzian_lines(FWHM, freq, res_frq) result(gamma)
          double precision, intent(in) :: res_frq, FWHM, freq
          double precision :: gamma
          gamma = ((FWHM)**2)/((FWHM)**2 + 4*(freq - res_frq)**2)
     end function lorentzian_lines

     function resonant_freq(q, t_rt) result(frq)
          double precision, intent(in) :: q, t_rt
          double precision :: frq
         frq = q/t_RT
     end function
         

!*****************************************************************************************
!create_spectrum will make a spectrum for feeding into the fabry-perot interferometer
!temp is temperature in kelvin
!peaks is a character array with information for the absorption peaks 
!        expects peaks in format:
!       peaks(1) = 'location = 352e-9 width = 10e-9 size = 0.8'
!       peaks(2) = 'location = 568e-9 width = 5e-9 size = 0.3'
!         location and width are in metres and size is from 0-1 and represents
!         percentage drop in spectrum. ie. a peak of 0.8 is 80% of the original height
!
!noise is how noisy the spectrum should be from 0 to 1. eg. noise = 0.01 is a spectrum with
!noise 1% of the height of the spectrum
!
!min, max are the bounds of the spectrum in m, eg. 300e-9 and 600e-9
!resolution is the step between each point eg. 1e-9
!
!******************************************************************************************

     function create_spectrum(temp, peaks, noise, min, max, resolution) result(spectrum)
          double precision, intent(in) :: temp, &
                                          noise, &
                                          min, &
                                          max, &
                                          resolution
          double precision, intent(in) :: peaks(:,:)
          double precision, dimension(:,:), allocatable :: spectrum
          double precision, dimension(3) :: peak_inf
          integer :: i, k, counter
          integer :: peak_num
          double precision :: frequency, z        

          counter = nint((max-min)/resolution)

          peak_num = size(peaks, 1)
          allocate (spectrum(2,counter)) 
          do k=1,counter
             z = (k-1)*resolution + min
             spectrum(1, k) = z
             frequency = c/z
             spectrum(2, k) = add_noise(noise, planck_spectrum(frequency, temp))
          end do
          do i=1,peak_num
             call add_peaks(spectrum, peaks(i,1), peaks(i,2), peaks(i,3), counter)
          end do
            
          end function create_spectrum
          
           function planck_spectrum(f, T) result(intensity)
               double precision, intent(in) :: f, T
               double precision :: intensity
               !intensity = 2*kb*T*((f**2)/(c**2))*(((h*f)/(kb*T))/(exp((h*f)/(kb*T))-1))
                intensity = ((2*h*f**3)/(c**2))*(1/(exp((h*f)/(kb*T))-1))
           end function planck_spectrum

           function add_noise(noise, val) result(new_val)
               double precision, intent(in) :: noise, val
               double precision :: new_val, noise_factor, noise_change
               real :: r1, r2
               call random_number(r1)
               call random_number(r2)
               noise_factor = noise*val
               noise_change = r1 * noise_factor
               if (r2 <= 0.5) then
                  noise_change = noise_change*(-1)
               end if            
               new_val = noise_change + val 
            end function add_noise
         
            subroutine add_peaks(spectrum, peak, width, peak_size, spectrum_size)
                integer, intent(in) :: spectrum_size 
                double precision, intent(inout) :: spectrum(:,:)
                double precision, intent(in) :: peak, width, peak_size
                integer :: peak_index, width_indicies, i
                double precision :: resol
                resol = abs(spectrum(1,2) - spectrum(1,1))
                width_indicies = nint(width/resol)
                peak_index = locate_peaks(spectrum(1,:), peak, size(spectrum, 2))
                do i=peak_index-width_indicies,peak_index+width_indicies
                    if (i==peak_index) then
                        spectrum(2,i) = spectrum(2,i) * peak_size
                    else if (i<peak_index) then
                         spectrum(2,i) = spectrum(2,i) * (1-(width_indicies-(peak_index-i))*((1-peak_size)/width_indicies))
                    else if (i>peak_index) then
                         spectrum(2,i) = spectrum(2,i) * (1-(width_indicies-(i-peak_index))*((1-peak_size)/width_indicies))
                    end if
                end do
    
            end subroutine add_peaks
              function locate_peaks(wavelengths, peak_loc, array_size) result(wavelength_loc)
                 integer, intent(in) :: array_size
                 double precision, dimension(array_size), intent(in) :: wavelengths
                 double precision, intent(in) :: peak_loc
                 integer, dimension(1) :: wavelength_loc1
                 integer :: wavelength_loc
                 double precision, dimension(array_size) :: wavelength2
                 wavelength2 = abs(wavelengths - peak_loc)
                 wavelength_loc1 =  findloc(wavelength2, minval(wavelength2))
                 wavelength_loc = wavelength_loc1(1)   
              end function locate_peaks  
           
             function airy_circ(phase, R1, R2) result(airy)
                  double precision, intent(in) :: R1, R2, phase
                  double precision :: airy
                  airy = 1/((1-sqrt(R1*R2))**2 + 4*sqrt(R1*R2)*sin(phase)**2)
             end function airy_circ
             
             function airy_other(phase, R1, R2, section) result(airy)
                  double precision, intent(in) :: R1, R2, phase
                  character (len = *) :: section
                  double precision :: airy, factor
                  if (trim(section) == 'b_circ') then
                     factor = R2
                  else if (trim(section) == 'trans') then
                     factor = 1. - R2
                  else if (trim(section) == 'back') then
                     factor = 1. - R1
                  else if (trim(section) == 'emit') then
                     factor = 1. - R1*R2
                  else if (trim(section) == 'circ_prime') then
                     factor = 1. - R1 
                  else if (trim(section) == 'b_circ_prime') then
                     factor = (1. - R1)*R2
                  else if (trim(section) == 'trans_prime') then
                     factor = (1. - R1)*(1. - R2)             
                  else if (trim(section) == 'back_prime') then
                     factor = (1. - R1)**2
                  else if (trim(section) == 'emit_prime') then
                     factor = (1. - R1)*(1.-R1*R2)
                  end if
                  airy = airy_circ(phase, R1, R2) * factor
               end function airy_other
            function gamma_trans_prime(freq, res_freq, length, R1, R2) result(gamma)
                  double precision, intent(in) :: freq, res_freq, length, R1, R2
                  double precision :: tau, fwhm, t_rt, gamma
                  t_rt = time_round_trip(length)
                  tau = photon_decay_time(R1, R2, t_rt)
                  fwhm = lorentzian_FWHM(tau)
                  gamma = (((1.-R1)*(1.-R2))/(1-R1*R2))*spectral_line_shape(res_freq, freq, fwhm)
             end function gamma_trans_prime                    

             function scan_cavity(spectrum_in, R1, R2, length_init, start_lamda, finish_lamda, step) result(spectrum_out)
!                   =====================================================================================================
                    double precision, intent(in) :: R1, R2, length_init, start_lamda, finish_lamda, step
                    double precision, intent(in) :: spectrum_in(:,:)
                    double precision, dimension(:,:), allocatable :: spectrum_out
                    double precision :: start_point, end_point, length_actual, current_length, current_trt, current_phase, &
                                        I_inc, freq, I_trans, total_trans = 0, bandwidth
                    integer :: counter, counter2, i, j
!                   ======================================================================================================
                    start_point = nint((length_init/start_lamda))*start_lamda
                    end_point = start_point + (finish_lamda - start_lamda)                    
                    counter = nint((end_point - start_point)/step)
                    counter2 = size(spectrum_in,2)
                    allocate (spectrum_out(3,counter))
                    bandwidth = spectrum_in(1,2) - spectrum_in(1,1)
                    do i=1,counter
                          total_trans = 0
                          current_length = start_point  + i*step
                          spectrum_out(1,i) = current_length
                          spectrum_out(2,i) = length2lamda(current_length, start_point, start_lamda)
                          current_trt = time_round_trip(current_length)
                          do j=1,counter2
                             freq =  c/spectrum_in(1,j)
                             I_inc =  spectrum_in(2,j)                           
                            ! current_phase = phase_shift(freq, current_trt)
                            ! I_trans = airy_other(current_phase, R1, R2, 'trans_prime')*I_inc*bandwidth
                             I_trans = gamma_trans_prime(freq, c/(start_lamda + i*step), current_length, R1, R2)*I_inc
                             total_trans = total_trans + I_trans
                             I_trans = 0
                          end do
                          spectrum_out(3,i) = total_trans * bandwidth
                     end do 
                end function scan_cavity
                
                function length2lamda(length, start_point, start_lamda) result(lamda)
                    double precision, intent(in) :: length, start_lamda, start_point
                    double precision :: q, fsr, t_rt, lamda
                    lamda = (length - start_point) + start_lamda
                end function length2lamda
            
end module fp_sim      
