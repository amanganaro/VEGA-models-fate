import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.persistence_soil_irfmn.ismPersistenceSoilIrfmn;
import model.ModelExecutionTest;

public class PersistenceSoilIRFMNTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismPersistenceSoilIrfmn();
    }
}
