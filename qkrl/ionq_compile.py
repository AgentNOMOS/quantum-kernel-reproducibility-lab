"""Local, offline lowering to the IonQ Forte-family native gate set {GPI, GPI2, ZZ} with the
gate classes and equivalence library shipped in qiskit-ionq. No IonQ API call is made: the
Target is built locally (qiskit-ionq's own backend object would query the API for the qubit
count). The result is what WE would submit; IonQ's server-side handling of non-native (QIS)
submissions is not observable offline.
"""
from qiskit import transpile
from qiskit.circuit import Parameter
from qiskit.circuit.library import Measure
from qiskit.transpiler import Target
from qiskit_ionq import GPIGate, GPI2Gate, ZZGate
from qiskit_ionq.ionq_equivalence_library import add_equivalences

add_equivalences()


def forte_native_target(num_qubits=9):
    """Mirror of qiskit_ionq IonQBackend._make_target() for gateset='native', family 'forte'."""
    tg = Target(num_qubits=num_qubits)
    tg.add_instruction(GPIGate(Parameter("phi")))
    tg.add_instruction(GPI2Gate(Parameter("phi")))
    tg.add_instruction(ZZGate(Parameter("theta")))
    tg.add_instruction(Measure())
    return tg


def to_forte_native(qc, optimization_level=3, seed=20261014):
    return transpile(qc, target=forte_native_target(qc.num_qubits),
                     optimization_level=optimization_level, seed_transpiler=seed)
