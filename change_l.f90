program change_l
    use fp_sim
    use rsplot
    use rsutil
    real :: reflect_1 = 0.95, &
            reflect_2 = 0.98
    double precision :: length=0.3, increment=5e-9, x, y
    double precision, dimension(:,:), allocatable :: spectrum_in, spectrum_out
    double precision :: temp=6000., noise=0.005, min=600e-9, max=1100e-9, resolution=0.1e-9
    double precision :: R1=0.95, R2=0.95, start_lamda = 600e-9, finish_lamda = 1050e-9, step = 1e-9
    double precision, dimension(6,3) :: peaks
    integer :: counter = 0, counter2
    real :: upper_in, lower_in, lower_out
    real, dimension(1000) :: upper_out, length_init
    real, dimension(:), allocatable :: spect_x_in, spect_y_in, spect_x_out, spect_y_out
    double precision :: end_point, start_point

!===========================================================================================
!    Create spectrum for input
!      
!===========================================================================================
     
     peaks(1,:) = (/627d-9, 2d-9, 1d-1/)
     counter = nint((max-min)/resolution)
     allocate (spectrum_in(2,counter))
     allocate (spect_x_in(counter))
     allocate (spect_y_in(counter))
     spectrum_in = create_spectrum(temp=temp, peaks=peaks, noise=noise, min=min, max=max, resolution=resolution)
     spect_x_in = real(spectrum_in(1,:))
     spect_y_in = real(spectrum_in(2,:))
     upper_in = maxval(spect_y_in)
     lower_in = minval(spect_y_in)
!===========================================================================================
!    Simulate Fabry-Perot Spectrometer
!
!===========================================================================================
    do i=1,1000
     length_init = vector(1000, 0.01, 0.3)
     start_point = nint((length_init(i)/start_lamda))*start_lamda
     end_point = start_point + (finish_lamda - start_lamda)                                                
     counter2 = nint((end_point - start_point)/step)
     allocate (spectrum_out(3,counter2))
     allocate (spect_x_out(counter2))
     allocate (spect_y_out(counter2))
     spectrum_out = scan_cavity(spectrum_in, R1, R1, length_init(i)*1D0, start_lamda, finish_lamda, step)
     spect_x_out = real(spectrum_out(2,:))
     spect_y_out = real(spectrum_out(3,:))
     upper_out(i) = maxval(spect_y_out)/upper_in
     lower_out = minval(spect_y_out)
     deallocate (spectrum_out)
     deallocate (spect_x_out)
     deallocate (spect_y_out)
    end do
!============================================================================================
!    Plot the two spectra
!
!============================================================================================    
     call ps_init(pform='A0_P', rot90=0, fileout='changing_l')
     call ps_frame(1, length_init(1), minval(upper_out), length_init(1000), maxval(upper_out), &
                   xbl=40.,ybl=650. , clip=1, &
                   fill=ps_color([.94,.94,.94]))
     call ps_axis(1, ax='Xx', glp=0, d_tm=5., title= &
          'Wavelength [m]')
     call ps_axis(1, ax='Yy', title='Intensity Reduction Factor', glp=0)
     call ps_plot(1, 1000, length_init, upper_out, ci=1)    
     call ps_exit    


end program change_l
