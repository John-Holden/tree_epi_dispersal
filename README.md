# Tree epidemiology: a variable dispersal model
<h3> A modelling approach to tree epidemiology: a wind-borne dispersal model. </h3>
<p>A three parameter model of disease spread: tree density (\rho), dipersal distance (\ell), dispersal shape (infectivty) (\beta).
This implementation is general with a variable dispersal kernel.</p>

<h3> Directories  </h3>
<ol>
  <li>anim_dat : dir, store all animation data, methods called by model in and after simulation.</li>
  <li>ensemble_dat: dir, to store ensemble data. Move to data_store after to save pushing extra information.</li>
  <li>SIR: dir, holds SIR fitting least-squares and different implementations of SIR equations.</li>
   <li>plots: dir, hold plotting methods. 
 </li>
</ol>

<h4>Scripts</h4>
<ol>
  <li> run_HPC.py : execute on arc super computer.</li>
  <li> run_model.py : high-level execution of the model, either single simulations or ensembles.</li>
  <li> run_model.py : high-level execution of the model, either single simulations or ensembles.</li>
  <li> model.py : contains the variable dispersal model. </li>
</ol>