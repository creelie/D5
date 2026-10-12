import D5Lean.Cert.All
import D5Lean.Sound.Facts

/-!
# The certificate theorem

Each of the 34 listed matrices `vtx i` satisfies `VFacts`, and each listed edge `e` at `vtx i`
satisfies `EFacts`; a bounded edge ends at one of the 34 matrices, up to `Γ`.
-/

namespace D5

theorem vertex_ok (i : Nat) (hi : i < 34) : (vtx i).ok = true := by
  have h := vertices_ok
  rw [List.all_eq_true] at h
  apply h
  unfold vtx
  rw [List.getD_eq_getElem _ _ (by rw [vertices_length]; exact hi)]
  exact List.getElem_mem _

/-- The finite part of the proof: the vertices and the edges of the certificate. -/
theorem certificate :
    (∀ i < 34, VFacts (vtx i)) ∧
    ∀ i < 34, ∀ e ∈ edgesOf i, (e.bounded = true → e.target < 34) ∧
      EFacts (vtx i) (vtx (if e.bounded then e.target else i)) e := by
  refine ⟨fun i hi => vertex_sound _ (vertex_ok i hi), fun i hi e he => ?_⟩
  have h := edges_ok i hi
  rw [List.all_eq_true] at h
  have h1 := h e he
  simp only [edgeOK, Bool.and_eq_true] at h1
  obtain ⟨h2, h3⟩ := h1
  have hs : shape (vtx i).J 6 = true := by
    have := vertex_ok i hi
    simp only [Vtx.ok, Bool.and_eq_true] at this
    exact this.1.1.1.1.1
  refine ⟨fun hb => ?_, edge_sound _ _ e hs h3⟩
  simpa [hb] using h2

end D5
